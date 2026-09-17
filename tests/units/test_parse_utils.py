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


def _run_argv(tmp_path, *extra):
    """Minimal valid argv for parse_args; every path it takes must exist."""
    anno = tmp_path / "in.anno"; anno.touch()
    bed = tmp_path / "target.bed"; bed.touch()
    return [
        "-i", anno.as_posix(),
        "-b", bed.as_posix(),
        "-o", tmp_path.as_posix(),
        "-n", "sample",
        *extra,
    ]


def test_off_mode_defaults_off(tmp_path):
    from SOPRANO.utils import parse_utils

    ns = parse_utils.parse_args(_run_argv(tmp_path))
    assert ns.off_mode is False


def test_off_mode_flag_sets_it(tmp_path):
    from SOPRANO.utils import parse_utils

    ns = parse_utils.parse_args(_run_argv(tmp_path, "--off_mode"))
    assert ns.off_mode is True


def test_off_mode_selects_the_min30_length_files(tmp_path):
    """The whole point of the flag reaching GlobalParameters."""
    from SOPRANO.core import objects
    from SOPRANO.utils import parse_utils

    plain = objects.TranscriptPaths.defaults()
    min30 = objects.TranscriptPaths.defaults(min30=True)

    # Separate directories: GlobalParameters caches its parameters and
    # refuses to reuse a directory whose cached run differs -- off_mode
    # included, which is itself worth knowing.
    on_dir = tmp_path / "on"; on_dir.mkdir()
    off_dir = tmp_path / "off"; off_dir.mkdir()

    ns_on = parse_utils.parse_args(_run_argv(on_dir, "--off_mode"))
    params_on = objects.GlobalParameters.from_namespace(ns_on)
    assert params_on.off_mode is True
    assert params_on.transcripts.transcript_length == min30.transcript_length
    assert (
        params_on.transcripts.protein_transcript_length
        == min30.protein_transcript_length
    )

    ns_off = parse_utils.parse_args(_run_argv(off_dir))
    params_off = objects.GlobalParameters.from_namespace(ns_off)
    assert params_off.off_mode is False
    assert params_off.transcripts.transcript_length == plain.transcript_length


def test_explicit_transcript_paths_survive_off_mode(tmp_path):
    """A path the user named is theirs, mode or no mode."""
    from SOPRANO.core import objects
    from SOPRANO.utils import parse_utils

    mine = tmp_path / "my_own.length"; mine.touch()
    ns = parse_utils.parse_args(
        _run_argv(tmp_path, "--off_mode", "-t", mine.as_posix())
    )
    params = objects.GlobalParameters.from_namespace(ns)
    assert params.transcripts.transcript_length == mine
