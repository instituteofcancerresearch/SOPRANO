import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity

from SOPRANO.utils.print_utils import task_output

_BANDWIDTH_POW_MIN = -3
_BANDWIDTH_POW_MAX = 2
_BANDWIDTH_N_DENSITY = 500


def sanitize_sklearn_input(
    data: pd.DataFrame, key: str, make_2d=True, sort=False
):
    numpy_array = data[key].to_numpy().copy()

    if sort:
        numpy_array.sort()

    if make_2d:
        return numpy_array[:, None]
    else:
        return numpy_array


def build_gaussian_kde(
    data: pd.DataFrame,
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

    input_values = sanitize_sklearn_input(data, key)

    grid.fit(input_values)

    task_output(f"Optimal bandwidth: {grid.best_params_}")

    return grid.best_estimator_


def estimate_density(data: pd.DataFrame, key: str, return_estimator=False):
    estimator = build_gaussian_kde(data, key)

    data_1d = np.linspace(0.25 * data[key].min(), data[key].max() * 1.75, 1000)
    data_2d = data_1d[:, None]

    lnp = estimator.score_samples(data_2d)

    p = np.exp(lnp)

    if return_estimator:
        return data_1d, p, estimator

    return data_1d, p
