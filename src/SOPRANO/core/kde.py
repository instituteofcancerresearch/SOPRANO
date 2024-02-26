import warnings
from typing import Callable

import numpy as np
import pandas as pd
from scipy.integrate import quad
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity

from SOPRANO.utils.print_utils import task_output

_BANDWIDTH_POW_MIN = -3
_BANDWIDTH_POW_MAX = 2
_BANDWIDTH_N_DENSITY = 500


def sanitize_sklearn_input(
    data_frame: pd.DataFrame, key: str, make_2d=True, sort=False
):
    numpy_array = data_frame[key].to_numpy().copy()

    if sort:
        numpy_array.sort()

    if make_2d:
        return numpy_array[:, None]
    else:
        return numpy_array


def build_gaussian_kde(
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

    input_values = sanitize_sklearn_input(null_hypothesis_samples, key)

    grid.fit(input_values)

    task_output(f"Optimal bandwidth: {grid.best_params_}")

    return grid.best_estimator_


def estimate_density(
    null_hypothesis_samples: pd.DataFrame, key: str, return_estimator=False
):
    estimator = build_gaussian_kde(null_hypothesis_samples, key)

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


def probability_estimator(
    null_hypothesis_samples: pd.DataFrame,
    key: str,
    pow_min=_BANDWIDTH_POW_MIN,
    pow_max=_BANDWIDTH_POW_MAX,
    n_density=_BANDWIDTH_N_DENSITY,
):
    base_estimator = build_gaussian_kde(
        null_hypothesis_samples, key, pow_min, pow_max, n_density
    )

    def _estimator_to_vectorize(x):
        return np.exp(base_estimator.score_samples(np.array([x])[:, None]))

    return np.vectorize(_estimator_to_vectorize)


def renormalize_estimator(
    null_hypothesis_estimator: Callable,
):
    warnings.warn(
        "Estimator requires renormalization! "
        "Currently unimplemented. "
        "pvalues will be unreliable."
    )
    return null_hypothesis_estimator


def determine_integration_bounds(
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

    def find_convergence(first_step, method: Callable, _iterations=0):
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
            return find_convergence(
                next_step, method, _iterations=_iterations + 1
            )

        return first_step

    lower_integral_bound = find_convergence(
        sample_values_min, lambda x: x - step_distance
    )

    upper_integral_bound = find_convergence(
        sample_values_min, lambda x: x + step_distance
    )

    return lower_integral_bound, upper_integral_bound


def estimate_pvalue(
    estimator: Callable,
    observed_value: float,
    lower_integral_bound: float,
    upper_integral_bound: float,
):
    pvalue_left = quad(estimator, lower_integral_bound, observed_value)[0]
    pvalue_right = quad(estimator, observed_value, upper_integral_bound)[0]

    return pvalue_left, pvalue_right


def estimate_std(
    sample_values: np.ndarray,
):
    return sample_values.std()
