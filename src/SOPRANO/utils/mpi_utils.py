try:
    from mpi4py import MPI
except ImportError:
    import _MPI as MPI

print(f"-- Loaded MPI module: {MPI.__file__}")
