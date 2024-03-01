import pathlib
import tempfile

from SOPRANO.core import analysis


def test__test_sum_possible_across_region(_sum_possible_across_region_fixture):
    (
        mock_sites_input,
        expected_sum_output,
    ) = _sum_possible_across_region_fixture

    with tempfile.TemporaryDirectory() as _tmpdir:
        _tmpdir = pathlib.Path(_tmpdir)

        temp_file = _tmpdir / "mock_data.sites"

        sum_file = pathlib.Path(_tmpdir) / "mock_data.sums"

        with open(temp_file, "w") as f:
            f.write(mock_sites_input)

        analysis._sum_possible_across_region(temp_file, sum_file)

        with open(sum_file, "r") as f:
            computed = f.readlines()

        computed = [line.strip() for line in computed]
        expected = [line.strip() for line in expected_sum_output.split("\n")]

        for _computed, _expected in zip(computed, expected):
            assert _computed == _expected
