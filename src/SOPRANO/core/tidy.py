import pathlib

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

    expression.append(r"\)")

    return expression
