try:
    from mpi4py import MPI
except ImportError:
    import _MPI as MPI

from copy import deepcopy
from typing import Callable

print(f"-- Loaded MPI module: {MPI.__file__}")


def item_selection(
    items_to_distribute: list | dict, mpi_rank: int, mpi_size: int
) -> list:
    if isinstance(items_to_distribute, dict):
        items_to_distribute = list(items_to_distribute.keys())

    if mpi_size == 1:
        return items_to_distribute

    items_to_distribute = deepcopy(items_to_distribute)

    current_idx = 0
    output: dict[int, list] = {k: [] for k in range(mpi_size)}

    while items_to_distribute:
        if current_idx == mpi_size:
            current_idx = 0

        output[current_idx].append(items_to_distribute.pop(0))

        current_idx += 1

    return output[mpi_rank]


def as_single_process(
    current_proc: int, mpi_comm: MPI.COMM_WORLD, executing_proc: int = 0
):
    def decorator(function: Callable):
        def wrapper(*args, **kwargs):
            if current_proc == executing_proc:
                result = function(*args, **kwargs)
            else:
                result = None
            mpi_comm.Barrier()
            return result

        return wrapper

    return decorator
