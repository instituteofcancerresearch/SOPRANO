import pathlib
import tempfile

import pytest

import SOPRANO.prepare_coordinates as prep_coords
from SOPRANO.objects import AnalysisPaths, TranscriptPaths


def tab_line(*args):
    return "\t".join([str(arg) for arg in args]) + "\n"


# Snippet from TCGA-05-4396-01A-21D-1855-08.annotated
mock_input_content = [
    tab_line(
        "10_118957090_G/C",
        "10:118957090",
        "C",
        "ENSG00000186795",
        "ENST00000334549",
        "Transcript",
        "missense_variant",
        "91",
        "91",
        "31",
        "V/L",
        "Gtg/Ctg",
        "-",
        "IMPACT=MODERATE;STRAND=1;SYMBOL=KCNK18;SYMBOL_SOURCE=HGNC;HGNC_ID=19439",
    ),
    tab_line(
        "10_16528530_G/A",
        "10:16528530",
        "A",
        "ENSG00000165983",
        "ENST00000378000",
        "Transcript",
        "synonymous_variant",
        "858",
        "612",
        "204",
        "R",
        "cgG/cgA",
        "-",
        "IMPACT=LOW;STRAND=1;SYMBOL=PTER;SYMBOL_SOURCE=HGNC;HGNC_ID=9590",
    ),
    tab_line(
        "10_22021940_G/T",
        "10:22021940",
        "T",
        "ENSG00000078403",
        "ENST00000307729",
        "Transcript",
        "missense_variant",
        "2509",
        "2331",
        "777",
        "Q/H",
        "caG/caT",
        "-",
        "IMPACT=MODERATE;STRAND=1;SYMBOL=MLLT10;SYMBOL_SOURCE=HGNC;HGNC_ID=16063",
    ),
    tab_line(
        "10_43015604_C/A",
        "10:43015604",
        "A",
        "ENSG00000234420",
        "ENST00000452075",
        "Transcript",
        "non_coding_transcript_exon_variant",
        "1991",
        "-",
        "-",
        "-",
        "-",
        "-",
        "IMPACT=MODIFIER;STRAND=-1;SYMBOL=ZNF37BP;SYMBOL_SOURCE=HGNC;HGNC_ID=13103",
    ),
    tab_line(
        "10_51363129_C/A",
        "10:51363129",
        "A",
        "ENSG00000225784",
        "ENST00000404618",
        "Transcript",
        "non_coding_transcript_exon_variant",
        "943",
        "-",
        "-",
        "-",
        "-",
        "-",
        "IMPACT=MODIFIER;STRAND=-1;SYMBOL=RP11-592B15.4;SYMBOL_SOURCE=Clone_based_vega_gene",
    ),
]

# Snippet from TCGA-05-4396.Expressed.IEDBpeps.SB.epitope.bed
mock_bed_content = [
    tab_line("ENST00000000233", 115, 124),
    tab_line("ENST00000000233", 164, 177),
    tab_line("ENST00000000412", 113, 124),
    tab_line("ENST00000001008", 27, 36),
    tab_line("ENST00000001008", 189, 198),
]

# Snippet from ensemble_transcript.length
mock_transcript_content = [
    tab_line("ENST00000000233", 543),
    tab_line("ENST00000000412", 834),
    tab_line("ENST00000000442", 1272),
    tab_line("ENST00000001008", 1380),
    tab_line("ENST00000001146", 1539),
]

# ensemble_transcript_protein.length
mock_protein_transcript_content = [
    tab_line("ENST00000000233", 180),
    tab_line("ENST00000000412", 277),
    tab_line("ENST00000000442", 423),
    tab_line("ENST00000001008", 459),
    tab_line("ENST00000001146", 512),
]


def check_expected_content(
    expected_content: list, written_content_path: pathlib.Path
):
    assert written_content_path.exists()

    with open(written_content_path, "r") as f:
        written_content = f.readlines()

    assert len(expected_content) == len(written_content)

    for e, w in zip(expected_content, written_content):
        assert e.strip() == w.strip()


@pytest.fixture
def test_files(tmp_path):
    inputs_dir: pathlib.Path = tmp_path.joinpath("inputs")
    transcripts_dir: pathlib.Path = tmp_path.joinpath("transcripts")
    tmpdir: pathlib.Path = tmp_path.joinpath("tmp")

    inputs_dir.mkdir(parents=True)
    transcripts_dir.mkdir(parents=True)
    tmpdir.mkdir(parents=True)

    anno_path = inputs_dir.joinpath("input.anno")
    bed_path = inputs_dir.joinpath("input.bed")

    trans_path = transcripts_dir.joinpath("transcript_length.txt")
    trans_prot_path = transcripts_dir.joinpath(
        "transcript_length_protein.length"
    )

    for _input_path, _input_content in zip(
        (anno_path, bed_path, trans_path, trans_prot_path),
        (
            mock_input_content,
            mock_bed_content,
            mock_transcript_content,
            mock_protein_transcript_content,
        ),
    ):
        with open(_input_path, "w") as f:
            f.writelines(_input_content)

    paths = AnalysisPaths("test_data", bed_path, tmpdir)
    transcripts = TranscriptPaths(trans_path, trans_prot_path)

    return paths, transcripts


@pytest.mark.dependency(name="_filter_transcript_file")
def test__filter_transcript_file(test_files):
    paths, transcripts = test_files

    expected_content = [
        tab_line("ENST00000000233", 543),
        tab_line("ENST00000000412", 834),
        tab_line("ENST00000001008", 1380),
    ]

    prep_coords._filter_transcript_file(
        paths.bed_path,
        transcripts.transcript_length,
        paths.filtered_transcript,
    )

    check_expected_content(expected_content, paths.filtered_transcript)


@pytest.mark.dependency(depends=["_filter_transcript_file"])
def test_filter_transcript_files(test_files):
    paths, transcripts = test_files

    expected_trans_content = [
        tab_line("ENST00000000233", 543),
        tab_line("ENST00000000412", 834),
        tab_line("ENST00000001008", 1380),
    ]

    expected_trans_protein_content = [
        tab_line("ENST00000000233", 180),
        tab_line("ENST00000000412", 277),
        tab_line("ENST00000001008", 459),
    ]

    prep_coords.filter_transcript_files(paths, transcripts)

    check_expected_content(expected_trans_content, paths.filtered_transcript)
    check_expected_content(
        expected_trans_protein_content, paths.filtered_protein_transcript
    )


def test__define_excluded_regions_for_randomization():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]
    expected_bed_content = [
        tab_line("chr1", 800, 1000),
        tab_line("chr1", 80, 180),
        tab_line("chr1", 1, 10),
        tab_line("chr1", 750, 10000),
        tab_line("chr1", 0, 2),
        tab_line("chr1", 0, 2),
        tab_line("chr1", 0, 2),
        tab_line("chr1", 0, 2),
    ]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)
        name = "test"
        bed_path = tmpdir.joinpath("bed.file")

        with open(bed_path, "w") as f:
            f.writelines(mock_bed_content)

        result_path = tmpdir.joinpath(f"{name}.exclusion.ori")

        paths = AnalysisPaths(name, bed_path, tmpdir)

        prep_coords._define_excluded_regions_for_randomization(paths)

        assert result_path.exists()

        with open(result_path, "r") as f:
            written_content = f.readlines()

        assert len(expected_bed_content) == len(written_content)

        for e, w in zip(expected_bed_content, written_content):
            assert e.strip() == w.strip()


def test__sort_excluded_regions_for_randomization():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]

    mock_filtered_protein = [
        tab_line("chr1", 10000),
        tab_line("chr2", 8000),
        tab_line("chr3", 5000),
        tab_line("chr2", 2000),
    ]

    input_exclusions = [
        tab_line("chr1", 300, 400),
        tab_line("chr1", 200, 250),
    ]

    # Default bedsort is by chrom then start pos
    expected_sorted_exclusions = [
        tab_line("chr1", 200, 250),
        tab_line("chr1", 300, 400),
    ]
    chr1_exclusions = [[200, 250], [300, 400]]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)
        name = "test"
        bed_path = tmpdir.joinpath("bed.file")
        paths = AnalysisPaths(name, bed_path, tmpdir)

        with open(paths.bed_path, "w") as f_bed:
            f_bed.writelines(mock_bed_content)

        with open(paths.exclusions, "w") as f_excl:
            f_excl.writelines(input_exclusions)

        with open(paths.filtered_protein_transcript, "w") as f_prot:
            f_prot.writelines(mock_filtered_protein)

        prep_coords._sort_excluded_regions_for_randomization(paths, seed=1234)

        assert paths.exclusions_sorted.exists()

        with open(paths.exclusions_sorted, "r") as f_sort:
            written_sorted_exclusions = f_sort.readlines()

        for e, w in zip(expected_sorted_exclusions, written_sorted_exclusions):
            assert e.strip() == w.strip()

        with open(paths.exclusions_shuffled, "r") as f_shuf:
            written_shuffled = f_shuf.readlines()

        for line in written_shuffled:
            line = line.strip()
            chrom, start, stop, *other = line.split("\t")
            if chrom == "chr1":
                start = int(start)
                stop = int(stop)
                for input_start, input_stop in chr1_exclusions:
                    overlap = (start <= input_stop) and (input_start <= stop)
                    assert not overlap


def test__randomize_with_target_file():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]

    mock_filtered_transcript = [
        tab_line("chr1", 10000),
        tab_line("chr2", 8000),
        tab_line("chr3", 5000),
        tab_line("chr2", 2000),
    ]

    mock_transcript_content = [
        tab_line("chr1", 2738),
        tab_line("chr2", 12),
    ]

    targets = [tab_line("chr1", 500, 1000), tab_line("chr1", 2000, 2500)]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)

        paths = AnalysisPaths(
            "test",
            tmpdir.joinpath("test.bed"),
            tmpdir,
            target_regions_path=tmpdir.joinpath("targets.bed"),
        )

        transcripts = TranscriptPaths(
            tmpdir.joinpath("trans.length"),
            tmpdir.joinpath("trans_prot.length"),
        )

        for path, lines in zip(
            (
                paths.bed_path,
                paths.filtered_transcript,
                paths.filtered_protein_transcript,
                paths.target_regions_path,
                transcripts.transcript_length,
                transcripts.protein_transcript_length,
            ),
            (
                mock_bed_content,
                mock_filtered_transcript,
                mock_filtered_transcript,
                targets,
                mock_transcript_content,
                mock_transcript_content,
            ),
        ):
            with open(path, "w") as f:
                f.writelines(lines)

        prep_coords._randomize_with_target_file(paths, transcripts, seed=1234)

        with open(paths.filtered_protein_transcript, "r") as f:
            written_filtered_lines = f.readlines()

        expected_filtered_lines = mock_filtered_transcript + [
            tab_line("chr1", 2738)
        ]

        for e, w in zip(expected_filtered_lines, written_filtered_lines):
            assert e.strip() == w.strip()


def test__non_randomized():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]

    expected_content = [
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 800, 1000, 24),
    ]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)

        paths = AnalysisPaths("test", tmpdir.joinpath("test.bed"), tmpdir)

        with open(paths.bed_path, "w") as f:
            f.writelines(mock_bed_content)

        prep_coords._non_randomized(paths)

        assert paths.exclusions_shuffled.exists()

        with open(paths.exclusions_shuffled, "r") as f:
            written = f.readlines()

        for e, w in zip(expected_content, written):
            assert e.strip() == w.strip()


def test__exclude_positively_selected_genes_disabled():
    pass
    # with tempfile.TemporaryDirectory() as _tmpdir:
    #     tmpdir = pathlib.Path(_tmpdir)
    #
    #     paths = AnalysisPaths("test", tmpdir.joinpath("test.bed"), tmpdir)
