from SOPRANO.core import objects


def test_Transcripts():
    defaults = objects.TranscriptPaths.defaults()

    assert defaults.transcript_length.exists()
    assert defaults.protein_transcript_length.exists()


def test_GRCh37():
    assert objects.GenomePaths.GRCh37().sizes.exists()


def test_GRCh38():
    assert objects.GenomePaths.GRCh38().sizes.exists()


def test_Transcripts_min30_paths():
    """OFF mode substitutes the >=30 amino acid length files.

    Only the two length files change; the transcript fasta is shared.
    Regenerate with scripts/make_min30_lengths.sh.
    """
    defaults = objects.TranscriptPaths.defaults()
    min30 = objects.TranscriptPaths.defaults(min30=True)

    assert min30.transcript_length.name == "ensemble_transcript_min30.length"
    assert (
        min30.protein_transcript_length.name
        == "ensemble_transcript_protein_min30.length"
    )
    assert min30.transcript_fasta == defaults.transcript_fasta
    assert min30.transcript_length != defaults.transcript_length

    # Now shipped, so they have to actually be there for OFF mode to run.
    assert min30.transcript_length.exists()
    assert min30.protein_transcript_length.exists()

    # One filter in two units: the protein file counts amino acids, the
    # transcript file bases, so both drop the same transcripts.
    n_protein = sum(1 for _ in min30.protein_transcript_length.open())
    n_transcript = sum(1 for _ in min30.transcript_length.open())
    assert n_protein == n_transcript
    assert n_protein < sum(1 for _ in defaults.protein_transcript_length.open())
