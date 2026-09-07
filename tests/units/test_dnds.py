import pandas as pd

import SOPRANO.core.dnds as dnds


def test__preprocess_dfs(tcga_05_4396_ssb192_cfg):
    (
        paths,
        merged_path,
        sites_extra_path,
        sites_intra_path,
    ) = tcga_05_4396_ssb192_cfg

    merged, sites_extra, sites_intra = dnds._preprocess_dfs(paths)

    assert merged.head().equals(pd.read_csv(merged_path, delimiter="\t"))
    assert sites_extra.head().equals(
        pd.read_csv(sites_extra_path, delimiter="\t")
    )
    assert sites_intra.head().equals(
        pd.read_csv(sites_intra_path, delimiter="\t")
    )


def test__compute_mutation_counts():
    mock_df = pd.DataFrame(
        {
            "EnsemblID": ["ENSTxxx"] * 10,
            "intronrate": [123] * 10,
            "extra_missense_variant": [0] * 10,
            "extra_synonymous_variant": [1] * 10,
            "intra_missense_variant": [0] * 10,
            "intra_synonymous_variant": [1] * 10,
        }
    )

    expected_series = pd.Series(
        {
            "extra_missense_variant": 0,
            "extra_synonymous_variant": 10,
            "intra_missense_variant": 0,
            "intra_synonymous_variant": 10,
            "mut_total_epitope": 10,
            "mut_total_non_epitope": 10,
        },
        index=[
            "extra_missense_variant",
            "extra_synonymous_variant",
            "intra_missense_variant",
            "intra_synonymous_variant",
            "mut_total_epitope",
            "mut_total_non_epitope",
        ],
    )

    computed_series = dnds._compute_mutation_counts(mock_df)

    assert computed_series.equals(expected_series)


def test__define_variables():
    mock_series = pd.Series({"m_1": 1, "m_2": 2})
    mock_df_1 = pd.DataFrame({"x_1": [3], "x_2": [4]})
    mock_df_2 = pd.DataFrame({"y_1": [5], "y_2": [6]})

    expected_series = pd.Series(
        {"m_1": 1, "m_2": 2, "x_1": 3, "x_2": 4, "y_1": 5, "y_2": 6}
    )

    computed_series = dnds._define_variables(mock_series, mock_df_1, mock_df_2)

    assert computed_series.equals(expected_series)


def test__rescale_intron_by_synonymous():
    mock_vars = pd.Series(
        {
            "mutsintron": 1,
            "intra_synonymous_variant": 1,
            "extra_synonymous_variant": 1,
            "intra_site_2": 1,
            "extra_site_2": 1,
        }
    )

    expected_value = 1.0
    assert dnds._rescale_intron_by_synonymous(mock_vars) == expected_value


def test__compute_kaks_intra_extra():
    mock_vars = pd.Series(
        {
            "intra_synonymous_variant": 1,
            "extra_synonymous_variant": 1,
            "intra_missense_variant": 2,
            "extra_missense_variant": 1,
            "intra_site_1": 1,
            "extra_site_1": 1,
            "intra_site_2": 2,
            "extra_site_2": 1,
        }
    )

    assert dnds._compute_kaks_extra(mock_vars) == 1
    assert dnds._compute_kaks_intra(mock_vars) == 4


def test__compute_kaks_intron():
    mock_vars = pd.Series(
        {
            "mutsintron": 1,
            "intra_synonymous_variant": 1,
            "extra_synonymous_variant": 1,
            "intra_missense_variant": 1,
            "extra_missense_variant": 1,
            "intra_site_1": 1,
            "extra_site_1": 1,
            "intra_site_2": 1,
            "extra_site_2": 1,
        }
    )

    assert dnds._compute_kaks_intron(mock_vars) == 1.0


def test__compute_conf_interval():
    # TODO: Implement
    pass


class _CoveragePaths:
    """Just enough of AnalysisPaths for _compute_coverage to run for real."""

    def __init__(self, tmp_path, off_mode, with_intron):
        self.off_mode = off_mode
        self.data_epitopes = tmp_path / "data_epitopes"
        self.epitope_nans = tmp_path / "epitope_nans"
        self.intra_epitope_nans = tmp_path / "intra_epitope_nans"
        self.intron_rate = tmp_path / "intron_rate"
        self.results_path = tmp_path / "results.tsv"

        # One transcript with both mutation classes, on and off target.
        self.data_epitopes.write_text(
            "ENST00000000233\t2\textra_missense_variant\n"
            "ENST00000000233\t4\textra_synonymous_variant\n"
            "ENST00000000233\t6\tintra_missense_variant\n"
            "ENST00000000233\t8\tintra_synonymous_variant\n"
        )
        self.epitope_nans.write_text("1000\t2000\n")
        self.intra_epitope_nans.write_text("3000\t4000\n")
        # An empty intron rate file means no intron correction, so no
        # Exonic_Intronic row -- which is the case OFF mode drops the p-value in.
        self.intron_rate.write_text(
            "ENST00000000233\t0.5\t3\t500\n" if with_intron else ""
        )


def _run_coverage(tmp_path, off_mode, with_intron):
    paths = _CoveragePaths(tmp_path, off_mode, with_intron)
    dnds._compute_coverage(paths)
    return pd.read_csv(paths.results_path, sep="\t", keep_default_na=False)


def test_off_mode_drops_pvalue_without_intron_correction(tmp_path):
    """calculateKaKsEpiCorrected_CI_mod4OFF.R prints the literal NA."""
    df = _run_coverage(tmp_path, off_mode=True, with_intron=False)
    assert list(df["Coverage"]) == ["Exonic_Only"]
    assert list(df["Pvalue"]) == ["NA"]


def test_off_mode_keeps_pvalue_with_intron_correction(tmp_path):
    """calculateKaKsEpiCorrected_CI_intron_V3_mod4OFF.R still prints Pval."""
    df = _run_coverage(tmp_path, off_mode=True, with_intron=True)
    assert list(df["Coverage"]) == ["Exonic_Only", "Exonic_Intronic"]
    assert "NA" not in list(df["Pvalue"])
    assert all(float(v) >= 0 for v in df["Pvalue"])


def test_pvalue_unaffected_outside_off_mode(tmp_path):
    for with_intron in (False, True):
        sub = tmp_path / f"intron_{with_intron}"
        sub.mkdir()
        df = _run_coverage(sub, off_mode=False, with_intron=with_intron)
        assert "NA" not in list(df["Pvalue"])
