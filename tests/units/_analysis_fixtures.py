import pytest

MOCK_SITES_INPUT = """ENST00000299004:11:348\tGAC_T\tGAC_T\t0\t1
ENST00000299004:1375:1825\tTAT_C\tTAT_C\t0\t0
ENST00000001657:696:991\tTGG_C\tTGG_C\t0\t0
ENST00000286386:923:1169\tAGC_A\tAGC_A\t0\t0
ENST00000079666:1281:1797\tACC_G\tACC_G\t0\t1
"""
EXPECTED_SUM_OUTPUT = """test_ACC_G_ACC_G\t0\t1
test_AGC_A_AGC_A\t0\t0
test_GAC_T_GAC_T\t0\t1
test_TAT_C_TAT_C\t0\t0
test_TGG_C_TGG_C\t0\t0
"""


@pytest.fixture
def _sum_possible_across_region_fixture():
    return MOCK_SITES_INPUT, EXPECTED_SUM_OUTPUT
