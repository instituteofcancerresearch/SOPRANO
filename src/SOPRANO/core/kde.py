import os
import warnings
from dataclasses import dataclass
from typing import Callable, Union

import numpy as np
import pandas as pd
from scipy.integrate import quad
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity

from SOPRANO.utils.print_utils import task_output

_BANDWIDTH_POW_MIN = -3
_BANDWIDTH_POW_MAX = 2
_BANDWIDTH_N_DENSITY = 500


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
    pow_min=_BANDWIDTH_POW_MIN,
    pow_max=_BANDWIDTH_POW_MAX,
    n_density=_BANDWIDTH_N_DENSITY,
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


def estimate_density(
    null_hypothesis_samples: pd.DataFrame, key: str, return_estimator=False
):
    estimator = _build_gaussian_kde(null_hypothesis_samples, key)

    data_1d = np.linspace(
        0.25 * null_hypothesis_samples[key].min(),
        null_hypothesis_samples[key].max() * 1.75,
        1000,
    )
    data_2d = data_1d[:, None]

    lnp = estimator.score_samples(data_2d)

    p = np.exp(lnp)

    if return_estimator:
        return data_1d, p, estimator

    return data_1d, p


def _probability_estimator(
    null_hypothesis_samples: pd.DataFrame,
    key: str,
    pow_min=_BANDWIDTH_POW_MIN,
    pow_max=_BANDWIDTH_POW_MAX,
    n_density=_BANDWIDTH_N_DENSITY,
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
    observed_value: float,
    lower_integral_bound: float,
    upper_integral_bound: float,
):
    pvalue_left = quad(estimator, lower_integral_bound, observed_value)[0]
    pvalue_right = quad(estimator, observed_value, upper_integral_bound)[0]

    return pvalue_left, pvalue_right


def _estimate_std_deviation(
    sample_values: np.ndarray,
):
    return sample_values.std()


def _data_available(null_hypothesis_samples_split: pd.DataFrame):
    return null_hypothesis_samples_split.shape[0] != 0


@dataclass
class EstimatorResults:
    estimator: Union[None, Callable]
    integration_bounds: Union[None, tuple]
    pvalues: Union[None, tuple]
    std: Union[None, float]

    @classmethod
    def null(cls):
        return cls(None, None, None, None)


_ABS_TOL = 1e-8
_REL_TOL = 1e-8
_STEP_PCT = 1e-2
_MAX_ITER = 1000

_ABS_TOL_ENV_VAR = "SOPRANO_KDE_ABS_TOL"
_REL_TOL_ENV_VAR = "SOPRANO_KDE_REL_TOL"
_STEP_PCT_ENV_VAR = "SOPRANO_KDE_STEP_PCT"
_MAX_ITER_ENV_VAR = "SOPRANO_KDE_MAX_ITER"


class _Data:
    def __init__(
        self,
        null_hypothesis_samples: pd.DataFrame,
        data: pd.DataFrame,
        column_key: str,
    ):
        self.null_hypothesis_samples = null_hypothesis_samples
        self.data_value = data[column_key].mean()
        self.column_key = column_key
        self.is_available = _data_available(null_hypothesis_samples)

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


class _HistogramData(_Data):
    def __init__(
        self,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
        column_key: str,
        sample_label: str,
        sample_color: str,
        sample_alpha: float,
        sample_hatch: str | None,
        sample_hist_type: str,
        data_label: str,
        data_color: str,
        data_line_style: str,
    ):
        super().__init__(
            null_hypothesis_samples=null_hypothesis_samples,
            data=data_results,
            column_key=column_key,
        )

        self.sample_histogram_kwargs = {
            "bins": "auto",
            "density": True,
            "label": sample_label,
            "facecolor": sample_color,
            "edgecolor": sample_color,
            "alpha": sample_alpha,
            "hatch": sample_hatch,
            "histtype": sample_hist_type,
        }

        self.sample_kde_kwargs = {
            "color": sample_color,
        }

        self.data_histogram_kwargs = {
            "label": data_label,
            "col": data_color,
            "ls": data_line_style,
        }

    @classmethod
    def on_exonic_only(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
        column_key: str,
        sample_label: str,
        sample_color: str,
        sample_alpha: float,
        sample_hatch: str | None,
        sample_hist_type: str,
        data_label: str,
        data_color: str,
        data_line_style: str,
    ):
        pass

    @classmethod
    def off_exonic_only(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
        column_key: str,
        sample_label: str,
        sample_color: str,
        sample_alpha: float,
        sample_hatch: str | None,
        sample_hist_type: str,
        data_label: str,
        data_color: str,
        data_line_style: str,
    ):
        pass

    @classmethod
    def on_exonic_intronic(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
        column_key: str,
        sample_label: str,
        sample_color: str,
        sample_alpha: float,
        sample_hatch: str | None,
        sample_hist_type: str,
        data_label: str,
        data_color: str,
        data_line_style: str,
    ):
        pass

    @classmethod
    def off_exonic_intronic(
        cls,
        null_hypothesis_samples: pd.DataFrame,
        data_results: pd.DataFrame,
        column_key: str,
        sample_label: str,
        sample_color: str,
        sample_alpha: float,
        sample_hatch: str | None,
        sample_hist_type: str,
        data_label: str,
        data_color: str,
        data_line_style: str,
    ):
        pass
