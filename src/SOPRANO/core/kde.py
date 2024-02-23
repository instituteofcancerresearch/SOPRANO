import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.neighbors import KernelDensity

_BANDWIDTH_POW_MIN = -2
_BANDWIDTH_POW_MAX = 2
_BANDWIDTH_N_DENSITY = 100


def build_gaussian_kde(
    data: pd.Series,
    pow_min=_BANDWIDTH_POW_MIN,
    pow_max=_BANDWIDTH_POW_MAX,
    n_density=_BANDWIDTH_N_DENSITY,
):
    bandwidths = 10 ** np.linspace(pow_min, pow_max, n_density)

    grid = GridSearchCV(
        KernelDensity(kernel="gaussian"),
        {"bandwidth": bandwidths},
        cv=LeaveOneOut(),
    )
    grid.fit(data.to_numpy()[:, None])

    gauss_kde = KernelDensity(
        kernel="gaussian", bandwidth=grid.best_params_["bandwidth"]
    )

    return gauss_kde
