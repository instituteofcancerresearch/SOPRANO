import pathlib
from argparse import Namespace
from unittest.mock import patch

from SOPRANO.utils.parse_utils import (
    fix_ns_species_arg,
    parse_args,
    parse_genome_args,
)


def test__fix_species_arg():
    ns = Namespace(species="Homo Sapiens")
    assert fix_ns_species_arg(ns).species == "homo_sapiens"


def test_parse_genome_args(capsys):
    with patch("sys.argv", ["parse_genome_args"]):
        args = parse_genome_args()
        assert args.species == "homo_sapiens"
        assert args.assembly == "GRCh38"
        assert args.release == 110
        assert args.primary_assembly is False
        assert args.download_only is False

    with patch(
        "sys.argv", ["parse_genome_args", "-s", "foo", "-a", "bar", "-p"]
    ):
        args = parse_genome_args()
        assert args.species == "foo"
        assert args.assembly == "bar"
        assert args.release == 110
        assert args.primary_assembly is True
        assert args.download_only is False


_NAME_FLAG = "-n"
_BED_FLAG = "-b"
_INPUT_FLAG = "-i"
_OUTPUT_FLAG = "-o"

_NAME_VALUE = pathlib.Path("dummy_name")
_BED_VALUE = pathlib.Path("dummy_bed")
_INPUT_VALUE = pathlib.Path("dummy_input")
_OUTPUT_VALUE = pathlib.Path("dummy_output")


def test_parse_args(capsys, tmp_path):
    TEST_ARGS = ["parse_args"]

    for flag, fname in zip(
        [_NAME_FLAG, _BED_FLAG, _INPUT_FLAG],
        [_NAME_VALUE, _BED_VALUE, _INPUT_VALUE],
    ):
        TEST_ARGS.append(flag)
        TEST_ARGS.append(tmp_path / fname)
        TEST_ARGS[-1].touch()
        TEST_ARGS[-1] = TEST_ARGS[-1].as_posix()

    for flag, dname in zip([_OUTPUT_FLAG], [_OUTPUT_VALUE]):
        TEST_ARGS.append(flag)
        TEST_ARGS.append(tmp_path / dname)
        TEST_ARGS[-1].mkdir()
        TEST_ARGS[-1] = TEST_ARGS[-1].as_posix()

    with patch("sys.argv", TEST_ARGS):
        args = parse_args()

        assert pathlib.Path(args.analysis_name).name == _NAME_VALUE.name
        assert pathlib.Path(args.bed_path).name == _BED_VALUE.name

        # Defaults... TODO: Improve!
        assert args.n_samples == 0
