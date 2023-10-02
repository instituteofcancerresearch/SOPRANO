import pandas as pd

from SOPRANO.core import objects
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.path_utils import Directories

input_file = Directories.examples("TCGA-05-4396-01A-21D-1855-08.annotated")
bed_file = Directories.immuno_humans(
    "TCGA-05-4396.Expressed.IEDBpeps.SB.epitope.bed"
)
name = "TCGA-05-4396"
genome_ref = "GRCh37"
exclude_drivers = True
release = 110

expected_path = Directories.int_tests("TCGA-05-4396.tsv")
assert expected_path.exists()


def test_pipeline(tmp_path):
    """
    Test the TCGA-05-4396 end-to-end for validation.

    :param tmp_path: invoked by pytest fixture - uses /tmp
    """
    params = objects.Parameters(
        analysis_name=name,
        input_path=input_file,
        bed_path=bed_file,
        cache_dir=tmp_path,
        random_regions=None,
        use_ssb192=False,
        use_random=False,
        exclude_drivers=exclude_drivers,
        seed=-1,
        transcripts=objects.TranscriptPaths.defaults(),
        genomes=objects.GenomePaths.GRCh37(),
    )

    run_pipeline(params)

    # Check filtered transcripts have been built
    assert params.filtered_transcript.exists()
    assert params.filtered_protein_transcript.exists()

    # Check drivers have (not) been excluded
    assert params.epitopes.exists()

    # Check epitope files been produced
    assert params.epitopes.exists()
    assert params.epitopes_cds.exists()

    # Check complement files been produced
    assert params.intra_epitopes_prot.exists()
    assert params.intra_epitopes_cds.exists()
    assert params.intra_epitopes_prot.exists()

    # Check that fasta files are built
    assert params.epitopes_cds_fasta.exists()
    assert params.intra_epitopes_cds_fasta.exists()

    assert params.epitopes_trans_regs.exists()
    assert params.intra_epitopes_trans_regs.exists()

    # Check theoretical site estimates computed
    assert params.epitopes_trans_regs.exists()
    assert params.intra_epitopes_trans_regs.exists()

    # Check that site summations have been performed
    assert params.epitopes_trans_regs_sum.exists()
    assert params.intra_epitopes_trans_regs_sum.exists()

    # Check col corrections have been applied to annotated input file
    assert params.sim_fixed.exists()
    assert params.col_corrected.exists()

    # Check contextualised calculation
    assert params.contextualised.exists()

    # Check flagged calculation
    assert params.flagged.exists()

    # Check triplet founts
    assert params.triplet_counts.exists()

    # Check that site corrections have been computed
    assert params.final_epitope_corrections.exists()
    assert params.final_intra_epitope_corrections.exists()

    # Check that intersect by frequency has been performed
    assert params.epitope_nans.exists()
    assert params.intra_epitope_nans.exists()

    # Check mutation counts
    assert params.variants_silent.exists()
    assert params.variants_nonsilent.exists()
    assert params.variants_missense.exists()
    assert params.variants_intronic.exists()
    assert params.raw_silent_count.exists()
    assert params.raw_nonsilent_count.exists()
    assert params.raw_missense_count.exists()
    assert params.in_silent_count.exists()
    assert params.in_nonsilent_count.exists()
    assert params.in_missense_count.exists()
    assert params.out_silent_count.exists()
    assert params.out_nonsilent_count.exists()
    assert params.out_missense_count.exists()

    # Check epitopes data file is build
    assert params.data_epitopes.exists()

    # Check intron rate has been computed
    assert params.intron_rate.exists()

    # Check results file exists!
    assert params.results_path.exists()

    computed_tsv = pd.read_csv(params.results_path, sep="\t")
    expected_tsv = pd.read_csv(expected_path, sep="\t")

    assert computed_tsv.equals(expected_tsv), (computed_tsv, expected_tsv)

    # coverage ON_dnds ON_lowci ON_highci ON_muts OFF_dnds OFF_lowci OFF_highci OFF_muts Pval ON_na ON_NA ON_ns ON_NS OFF_na OFF_NA OFF_ns OFF_NS                                                 # noqa: E501
    # ExonicOnly 0.170545315483698 0.0312367028034305 0.931140579583117 6 0.890687718057257 0.510646130660438 1.5535705238312 63 0.330510882590904 2 1974270 4 673405 46 19525700 17 6427220      # noqa: E501
    # ExonicIntronic 0.170545315483698 0.0312367028034305 0.931140579583117 6 0.890687718057257 0.510646130660438 1.5535705238312 63 0.330510882590904 2 1974270 4 673405 46 19525700 17 6427220  # noqa: E501