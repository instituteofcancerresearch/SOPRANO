import pathlib
import tempfile

from SOPRANO.core import analysis

mock_sites = """ENST00000299004:11:348\tGAC_T\tGAC_T\t0\t1
ENST00000299004:1375:1825\tTAT_C\tTAT_C\t0\t0
ENST00000001657:696:991\tTGG_C\tTGG_C\t0\t0
ENST00000286386:923:1169\tAGC_A\tAGC_A\t0\t0
ENST00000079666:1281:1797\tACC_G\tACC_G\t0\t1
"""

output = """test_ACC_G_ACC_G\t0\t1
test_AGC_A_AGC_A\t0\t0
test_GAC_T_GAC_T\t0\t1
test_TAT_C_TAT_C\t0\t0
test_TGG_C_TGG_C\t0\t0
"""


def test__test_sum_possible_across_region():
    with tempfile.TemporaryDirectory() as _tmpdir:
        _tmpdir = pathlib.Path(_tmpdir)

        temp_file = _tmpdir / "mock_data.sites"

        sum_file = pathlib.Path(_tmpdir) / "mock_data.sums"

        with open(temp_file, "w") as f:
            f.write(mock_sites)

        analysis._sum_possible_across_region(temp_file, sum_file)

        with open(sum_file, "r") as f:
            computed = f.readlines()

        computed = [line.strip() for line in computed]
        expected = [line.strip() for line in output.split("\n")]

        for _computed, _expected in zip(computed, expected):
            assert _computed == _expected
