import os
import pathlib
import warnings
from dataclasses import dataclass
from typing import Callable, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity

from SOPRANO.utils.print_utils import task_output

# Default values for gridsearch optimization & integration tols

_CV_LOG_MIN = -3
_CV_LOG_MAX = 2
_CV_LOG_STEPS = 250

_ABS_TOL = 1e-8
_REL_TOL = 1e-8
_STEP_PCT = 1e-2
_MAX_ITER = 1000

# Override by exporting the following variables

_CV_LOG_MIN_ENV_VAR = "SOPRANO_CV_LOG_MIN"
_CV_LOG_MAX_ENV_VAR = "SOPRANO_CV_LOG_MAX"
_CV_LOG_STEPS_ENV_VAR = "SOPRANO_CV_LOG_STEPS"

_ABS_TOL_ENV_VAR = "SOPRANO_KDE_ABS_TOL"
_REL_TOL_ENV_VAR = "SOPRANO_KDE_REL_TOL"
_STEP_PCT_ENV_VAR = "SOPRANO_KDE_STEP_PCT"
_MAX_ITER_ENV_VAR = "SOPRANO_KDE_MAX_ITER"


def _sanitize_sklearn_input(
    data_frame: pd.DataFrame, key: str, make_2d=True, sort=False
):
    numpy_array = data_frame[key].to_numpy().copy()

    if sort:
        numpy_array.sort()

    if make_2d:
        return numpy_array[:, None]
    else:
        return numpy_array


def _build_gaussian_kde(
    null_hypothesis_samples: pd.DataFrame,
    key: str,
    pow_min=float(os.getenv(_CV_LOG_MIN_ENV_VAR, _CV_LOG_MIN)),
    pow_max=float(os.getenv(_CV_LOG_MAX_ENV_VAR, _CV_LOG_MAX)),
    n_density=int(os.getenv(_CV_LOG_STEPS_ENV_VAR, _CV_LOG_STEPS)),
):
    bandwidths = 10 ** np.linspace(pow_min, pow_max, n_density)

    task_output("Performing grid search for KDE estimate")

    grid = GridSearchCV(
        KernelDensity(kernel="gaussian"),
        {"bandwidth": bandwidths},
        cv=LeaveOneOut(),
    )

    input_values = _sanitize_sklearn_input(null_hypothesis_samples, key)

    grid.fit(input_values)

    task_output(f"Optimal bandwidth: {grid.best_params_}")

    return grid.best_estimator_


def _probability_estimator(
    null_hypothesis_samples: pd.DataFrame,
    key: str,
    pow_min=_CV_LOG_MIN,
    pow_max=_CV_LOG_MAX,
    n_density=_CV_LOG_STEPS,
):
    base_estimator = _build_gaussian_kde(
        null_hypothesis_samples, key, pow_min, pow_max, n_density
    )

    def _estimator_to_vectorize(x):
        return np.exp(base_estimator.score_samples(np.array([x])[:, None]))

    return np.vectorize(_estimator_to_vectorize)


def _renormalize_estimator(
    null_hypothesis_estimator: Callable,
):
    warnings.warn(
        "Estimator requires renormalization! "
        "Currently unimplemented. "
        "pvalues will be unreliable."
    )
    return null_hypothesis_estimator


def _iterate_until_convergence(
    max_iterations,
    estimator,
    absolute_tol,
    relative_tol,
    first_step,
    method: Callable,
    _iterations=0,
):
    if _iterations > max_iterations:
        warnings.warn(
            f"Maximum number of iterations exceeded: {max_iterations}\n"
            f"exiting with P(dN/dS)={estimator(first_step)}"
        )

    next_step = method(first_step)

    first_prob = estimator(first_step)
    next_prob = estimator(next_step)

    absolute_diff = abs(first_prob - next_prob)
    relative_diff = absolute_diff / next_prob

    if absolute_diff > absolute_tol and relative_diff > relative_tol:
        return _iterate_until_convergence(
            max_iterations,
            estimator,
            absolute_tol,
            relative_tol,
            next_step,
            method,
            _iterations=_iterations + 1,
        )

    return first_step


def _determine_integration_bounds(
    estimator: Callable,
    sample_values_min: float,
    sample_values_max: float,
    absolute_tol: float = 1e-8,
    relative_tol: float = 1e-8,
    step_percent: float = 1e-2,
    max_iterations: int = 1000,
):
    samples_spread = sample_values_max - sample_values_min

    if sample_values_max < sample_values_min:
        raise ValueError(f"{sample_values_max} < {sample_values_min}")

    step_distance = step_percent * samples_spread

    lower_integral_bound = _iterate_until_convergence(
        max_iterations,
        estimator,
        absolute_tol,
        relative_tol,
        sample_values_min,
        lambda x: x - step_distance,
    )

    upper_integral_bound = _iterate_until_convergence(
        max_iterations,
        estimator,
        absolute_tol,
        relative_tol,
        sample_values_min,
        lambda x: x + step_distance,
    )

    return lower_integral_bound, upper_integral_bound


def _estimate_pvalue(
    estimator: Callable,
    observed_value: float | None,
    lower_integral_bound: float,
    upper_integral_bound: float,
) -> tuple[float, float]:
    pvalue_left = quad(estimator, lower_integral_bound, observed_value)[0]
    pvalue_right = quad(estimator, observed_value, upper_integral_bound)[0]

    return pvalue_left, pvalue_right


def _estimate_std_deviation(
    sample_values: np.ndarray,
):
    return sample_values.std()


def _samples_and_data_ara_available(
    null_hypothesis_samples_split: pd.DataFrame, data: pd.DataFrame
):
    data_available = True

    if null_hypothesis_samples_split.empty:
        warnings.warn("Null hypothesis samples are empty!")
        data_available = False

    if data.empty:
        warnings.warn("Data is empty!")
        data_available = False

    # return null_hypothesis_samples_split.shape[0] != 0 and not data.empty

    return data_available


def _null_estimator(*args, **kwargs):
    pass


@dataclass
class EstimatorResults:
    estimator: Callable
    integration_bounds: Union[None, tuple]
    pvalues: Union[None, tuple]
    std: Union[None, float]

    @classmethod
    def null(cls):
        return cls(_null_estimator, None, None, None)


class _Data:
    def __init__(
        self,
        null_hypothesis_samples: pd.DataFrame,
        data: pd.DataFrame,
        column_key: str,
    ):
        self.null_hypothesis_samples = null_hypothesis_samples
        self.data_value = None if data.empty else data[column_key].mean()

        # print(column_key)
        # print(self.data_value)
        # print(data)

        self.column_key = column_key
        self.is_available = _samples_and_data_ara_available(
            null_hypothesis_samples, data
        )

        if self.is_available:
            self.estimates = self._apply_estimator()
        else:
            self.estimates = EstimatorResults.null()

    def _apply_estimator(self) -> EstimatorResults:
        estimator = _probability_estimator(
            self.null_hypothesis_samples, self.column_key
        )

        samples = self.null_hypothesis_samples[self.column_key]

        integration_bounds = _determine_integration_bounds(
            estimator=estimator,
            sample_values_min=samples.min(),
            sample_values_max=samples.max(),
            absolute_tol=float(os.getenv(_ABS_TOL_ENV_VAR, _ABS_TOL)),
            relative_tol=float(os.getenv(_REL_TOL_ENV_VAR, _REL_TOL)),
            step_percent=float(os.getenv(_STEP_PCT_ENV_VAR, _STEP_PCT)),
            max_iterations=int(os.getenv(_MAX_ITER_ENV_VAR, _MAX_ITER)),
        )

        pvalues = _estimate_pvalue(
            estimator, self.data_value, *integration_bounds
        )

        std = _estimate_std_deviation(samples)

        return EstimatorResults(
            estimator=estimator,
            integration_bounds=integration_bounds,
            pvalues=pvalues,
            std=std,
        )


_EXONIC_BASE_PLOT_KWARGS = {
    "sample_label": "Exonic",
    "sample_color": "red",
    "sample_alpha": 0.33,
    "sample_hatch": None,
    "sample_hist_type": "stepfilled",
    "data_label": "Exonic Data",
    "data_color": "k",
    "data_line_style": "--",
}

_EXONIC_INTRONIC_BASE_PLOT_KWARGS = {
    "sample_label": "Exonic Intronic",
    "sample_color": "blue",
    "sample_alpha": 0.33,
    "sample_hatch": None,  # "/",
    "sample_hist_type": "stepfilled",  # "step",
    "data_label": "Exonic Intronic Data",
    "data_color": "tab:gray",
    "data_line_style": ":",
}


class _HistogramData(_Data):
    def __init__(
        self,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
        column_key: str,
        exonic_only: bool,
    ):
        if exonic_only:
            base_kwargs = _EXONIC_BASE_PLOT_KWARGS
        else:
            base_kwargs = _EXONIC_INTRONIC_BASE_PLOT_KWARGS

        super().__init__(
            null_hypothesis_samples=null_hypothesis_samples,
            data=data_results,
            column_key=column_key,
        )

        self.sample_hist_kwargs = {
            "bins": "auto",
            "density": True,
            "label": base_kwargs["sample_label"],
            "facecolor": base_kwargs["sample_color"],
            "edgecolor": base_kwargs["sample_color"],
            "alpha": base_kwargs["sample_alpha"],
            "hatch": base_kwargs["sample_hatch"],
            "histtype": base_kwargs["sample_hist_type"],
        }

        self.sample_plot_kwargs = {
            "color": base_kwargs["sample_color"],
        }

        self.data_vline_kwargs = {
            "label": base_kwargs["data_label"],
            "color": base_kwargs["data_color"],
            "ls": base_kwargs["data_line_style"],
        }

    @classmethod
    def on_exonic_only(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
    ):
        return cls(
            null_hypothesis_samples=null_hypothesis_samples,
            data_results=data_results,
            column_key="ON_dNdS",
            exonic_only=True,
        )

    @classmethod
    def off_exonic_only(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
    ):
        return cls(
            null_hypothesis_samples=null_hypothesis_samples,
            data_results=data_results,
            column_key="OFF_dNdS",
            exonic_only=True,
        )

    @classmethod
    def on_exonic_intronic(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
    ):
        return cls(
            null_hypothesis_samples=null_hypothesis_samples,
            data_results=data_results,
            column_key="ON_dNdS",
            exonic_only=False,
        )

    @classmethod
    def off_exonic_intronic(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
    ):
        return cls(
            null_hypothesis_samples=null_hypothesis_samples,
            data_results=data_results,
            column_key="OFF_dNdS",
            exonic_only=False,
        )

    def plot_hist(self, ax, zorder: int):
        if self.is_available:
            ax.hist(
                self.null_hypothesis_samples[self.column_key],
                zorder=zorder,
                **self.sample_hist_kwargs,
            )

    def plot_kde(self, ax, zorder: int):
        if self.is_available:
            x_space = np.linspace(*self.estimates.integration_bounds, 1000)
            ax.plot(
                x_space,
                self.estimates.estimator(x_space),
                zorder=zorder,
                **self.sample_plot_kwargs,
            )

    def plot_data(self, ax, zorder: int):
        if self.is_available:
            ax.axvline(
                self.data_value, zorder=zorder, **self.data_vline_kwargs
            )


class PlotData:
    def __init__(
        self,
        samples_exonic: pd.DataFrame,
        samples_exonic_intronic: pd.DataFrame,
        data_exonic: pd.DataFrame,
        data_exonic_intronic: pd.DataFrame,
    ):
        on_exonic = _HistogramData.on_exonic_only(samples_exonic, data_exonic)
        on_exonic_intronic = _HistogramData.on_exonic_intronic(
            samples_exonic_intronic, data_exonic_intronic
        )

        off_exonic = _HistogramData.off_exonic_only(
            samples_exonic, data_exonic
        )

        off_exonic_intronic = _HistogramData.off_exonic_intronic(
            samples_exonic_intronic, data_exonic_intronic
        )

        self.on_samples = [
            samples
            for samples in (on_exonic, on_exonic_intronic)
            if samples.is_available
        ]

        self.off_samples = [
            data
            for data in (off_exonic, off_exonic_intronic)
            if data.is_available
        ]

    def make_figure(self, job_cache: pathlib.Path):
        fig, (ax_on, ax_off) = plt.subplots(1, 2, figsize=(8, 4))

        ax_on.set_title("ON")
        ax_off.set_title("OFF")

        def update_value(old_value, new_value, method: Callable):
            return (
                new_value
                if old_value is None
                else method(old_value, new_value)
            )

        def render(data_sets, ax):
            x_min = None
            x_max = None

            for idx, data_set in enumerate(data_sets):
                # apply z order to retain sensible overlaying order in loop

                # data_set: _HistogramData

                data_set.plot_hist(ax=ax, zorder=idx)
                data_set.plot_kde(ax=ax, zorder=idx + 2)
                data_set.plot_data(ax=ax, zorder=idx + 4)

                x_min = update_value(
                    x_min, data_set.estimates.integration_bounds[0], min
                )
                x_max = update_value(
                    x_max, data_set.estimates.integration_bounds[1], max
                )

            ax.set_xlim(x_min, x_max)

        render(self.on_samples, ax_on)
        render(self.off_samples, ax_off)

        for ax in (ax_on, ax_off):
            ax.grid()
            ax.set_xlabel("$dN/dS$")
            ax.set_ylabel("$P(dN/dS)$")

        plt.tight_layout()

        ax_on.legend(
            frameon=False,
            bbox_to_anchor=(2 if len(self.on_samples) == 2 else 1.5, -0.2),
            ncol=4,
        )

        plt.savefig(job_cache.joinpath("figure.pdf"), bbox_inches="tight")
