import glob
import os
import tarfile

from SOPRANO.core.objects import Parameters

_PATTERNS_TO_IGNORE = ("*results.tsv", "*.log", "*.complete")


def find_files_to_tar_gz(params: Parameters):
    _cache_dir = params.cache_dir.as_posix()

    files_to_ignore = []

    for pattern in _PATTERNS_TO_IGNORE:
        files_to_ignore += glob.glob(f"{_cache_dir}/{pattern}")

    files_to_tar_gz = glob.glob(f"{_cache_dir}/*")

    # Remove ignored file pattens
    for _ in files_to_ignore:
        files_to_tar_gz.remove(_)

    return files_to_tar_gz


def tar_and_compress(params: Parameters):
    files_to_tar_gz = find_files_to_tar_gz(params)

    _tar_gz_path = params.cache_dir.joinpath("intermediate.data.tar.gz")

    if _tar_gz_path.exists():
        raise FileExistsError(_tar_gz_path)

    tar_gz_path = _tar_gz_path.as_posix()

    with tarfile.open(tar_gz_path, "w:gz") as tar:
        for path in files_to_tar_gz:
            print(f"adding {path} to tar archive")
            tar.add(path)

    for path in files_to_tar_gz:
        os.remove(path)
