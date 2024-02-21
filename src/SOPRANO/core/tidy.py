import pathlib

from SOPRANO.core.objects import Parameters
from SOPRANO.utils.sh_utils import pipe

_PATTERNS_TO_IGNORE = ("*results.tsv", "*.log")


def build_find_expression(
    cache_dir: pathlib.Path, patterns_to_ignore=_PATTERNS_TO_IGNORE
):
    expression = ["find", cache_dir.as_posix(), "-type", "f", "!", r"\("]

    list_of_patterns = list(patterns_to_ignore)

    n_patterns = len(list_of_patterns)

    while n_patterns > 0:
        pattern = list_of_patterns.pop(0)
        expression += ["-name", pattern]

        n_patterns = len(list_of_patterns)

        if n_patterns > 0:
            expression.append("-o")

    expression += [r"\)", "-print0"]

    return expression


def build_tar_gz_expression(cache_dir: pathlib.Path):
    tar_gz_path = cache_dir / "intermediate.data.tar.gz"

    return ["tar", "-cvzf", tar_gz_path.as_posix(), "--null", "-T", "-"]


def tar_and_compress(params: Parameters):
    _cache_dir = params.cache_dir.as_posix()

    find_expression = build_find_expression(params.cache_dir)
    tar_expression = build_tar_gz_expression(params.cache_dir)

    pipe(find_expression, tar_expression)
