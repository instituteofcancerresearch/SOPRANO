import pathlib
from contextlib import contextmanager, redirect_stdout
from io import StringIO
from time import time

import pandas as pd
import streamlit as st

from SOPRANO.core.objects import EnsemblData, Parameters
from SOPRANO.hla2ip import immunopeptidome_from_hla
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.parse_utils import fix_species_arg
from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe
from SOPRANO.utils.vep_utils import (
    _get_src_dst_link_pairs,
    _link_src_dst_pairs,
)


@contextmanager
def st_capture(output_func):
    """Taken from https://discuss.streamlit.io/t/
        cannot-print-the-terminal-output-in-streamlit/6602/2

    :param output_func: an instance of st.empty()

    """
    with StringIO() as stdout, redirect_stdout(stdout):
        old_write = stdout.write

        def new_write(string):
            ret = old_write(string)
            output_func(stdout.getvalue())
            return ret

        stdout.write = new_write
        yield


def _select_from_dict(selection: str, selection_dict: dict):
    selection_value = selection_dict[selection]
    st.text(f"Selected: {selection_value}")
    return selection_value


class _PipelineUI:
    @staticmethod
    def genome_reference(*args, **kwargs):
        pass

    @staticmethod
    def annotated_mutations(*args, **kwargs):
        pass

    @staticmethod
    def immunopeptidome(*args, **kwargs):
        pass

    @staticmethod
    def substitution_method(*args, **kwargs):
        pass

    @staticmethod
    def coordinates(*args, **kwargs):
        pass


class PipelineUIOptions(_PipelineUI):
    @staticmethod
    def genome_reference():
        homo_sapiens_dir = Directories.genomes_homo_sapiens()

        genome_dirs = [
            item for item in homo_sapiens_dir.glob("*") if item.is_dir()
        ]

        # Remove bad options (i.e. no toplevel fa and chrom files)
        for item in genome_dirs[::-1]:
            toplevel_path = item.glob("*dna*toplevel*.fa")
            chrom_path = item.glob("*dna*toplevel*.chrom")

            if len(list(toplevel_path)) == len(list(chrom_path)) == 1:
                pass
            else:
                genome_dirs.remove(item)

        genome_ids = [
            "{} - Ensembl release {}".format(*x.name.split("_")[::-1])
            for x in genome_dirs
        ]

        options_dict = {
            name: dir_path for name, dir_path in zip(genome_ids, genome_dirs)
        }

        return options_dict

    @staticmethod
    def annotated_mutations():
        options_dict = {}
        for directory in (
            Directories.examples(),
            Directories.app_annotated_inputs(),
        ):
            for x in directory.glob("*.anno*"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def immunopeptidome():
        options_dict = {}
        for directory in (
            Directories.immunopeptidomes_humans(),
            Directories.app_immunopeptidomes(),
        ):
            for x in directory.glob("*.bed"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def substitution_method():
        return {"SSB192": 192, "SSB7": 7}

    @staticmethod
    def coordinates():
        options_dict = {None: None}
        for x in Directories.app_coordinate_files().glob("*.bed"):
            options_dict[x.name] = x
        return options_dict


class PipelineUIProcessing(_PipelineUI):
    @staticmethod
    def genome_reference(genome_selection: str):
        assembly, release = genome_selection.split(" - Ensembl release ")
        data = EnsemblData(species="homo_sapiens", assembly=assembly)
        fasta_path = data.toplevel_fa_path(int(release))
        chrom_path = data.toplevel_chrom_path(int(release))
        st.text(f"Selected: {fasta_path}, {chrom_path}")
        return data.get_genome_reference_paths(int(release))

    @staticmethod
    def annotated_mutations(annotation_selection: str):
        options_dict = PipelineUIOptions.annotated_mutations()
        return _select_from_dict(annotation_selection, options_dict)

    @staticmethod
    def immunopeptidome(immunopeptidome_selection: str):
        options_dict = PipelineUIOptions.immunopeptidome()
        return _select_from_dict(immunopeptidome_selection, options_dict)

    @staticmethod
    def substitution_method(subs_selection: str):
        options_dict = PipelineUIOptions.substitution_method()
        return _select_from_dict(subs_selection, options_dict)

    @staticmethod
    def coordinates(coordinates_selection: str):
        options_dict = PipelineUIOptions.coordinates()
        return _select_from_dict(coordinates_selection, options_dict)

    @staticmethod
    def job_name(job_name: str):
        cache_dir = Directories.cache(job_name)
        st.text(f"Selected: {cache_dir}")
        return cache_dir


class _LinkVEPUI:
    @staticmethod
    def cache_location(*args, **kwargs):
        pass


class LinkVEPUIOptions(_LinkVEPUI):
    pass


class LinkVEPUIProcessing(_LinkVEPUI):
    @staticmethod
    def cache_location(cache_location: str):
        output = pathlib.Path(cache_location)
        st.text(f"Selected: {output}")
        return output


class _DownloaderUI:
    @staticmethod
    def species(*args, **kwargs):
        pass

    @staticmethod
    def assembly(*args, **kwargs):
        pass

    @staticmethod
    def release(*args, **kwargs):
        pass

    @staticmethod
    def type(*args, **kwargs):
        pass


class DownloaderUIOptions(_DownloaderUI):
    @staticmethod
    def type():
        return "toplevel", "primary_assembly"


class DownloaderUIProcessing(_DownloaderUI):
    @staticmethod
    def species(species_selection: str):
        output = fix_species_arg(species_selection)
        st.text(f"Selected: {output}")
        return output

    @staticmethod
    def assembly(assembly_selection: str):
        st.text(f"Selected: {assembly_selection}")
        return assembly_selection

    @staticmethod
    def release(release: str):
        output = int(release)
        st.text(f"Selected: {output}")

        if output > 110:
            st.text("[Warning] Oct 1 2023: Latest Ensembl release is 110")

        return output

    @staticmethod
    def type(type_selection: str):
        if type_selection not in ("toplevel", "primary_assembly"):
            raise ValueError(type_selection)

        st.text(f"Selected: {type_selection}")
        return type_selection


class _AnnotatorUI:
    pass


class AnnotatorUIOptions(_AnnotatorUI):
    pass


class AnnotatorUIProcessing(_AnnotatorUI):
    pass


class _ImmunopeptidomeUI:
    @staticmethod
    def hla_alleles(*args, **kwargs):
        pass

    @staticmethod
    def transcript_ids(*args, **kwargs):
        pass

    @staticmethod
    def subset_method(*args, **kwargs):
        pass

    @staticmethod
    def name(*args, **kwargs):
        pass


class ImmunopeptidomesUIOptions(_ImmunopeptidomeUI):
    @staticmethod
    def hla_alleles():
        hla_types_path = Directories.examples("TCGA_hlaTypesAll.tsv")
        options = pipe(
            ["cut", "-f3", hla_types_path.as_posix()],
            ["tr", ",", "\n"],
            ["sort", "-u"],
        ).split("\n")
        return options

    @staticmethod
    def transcript_ids():
        hla_binders_path = Directories.data(
            "allhlaBinders_exprmean1.IEDBpeps.bed.unique_ids"
        )

        with open(hla_binders_path, "r") as f:
            transcript_options = f.read()
        # eol marker generates empty so excluded
        return transcript_options.split("\n")[:-1]

    @staticmethod
    def subset_method():
        return "None", "Retention", "Exclusion"


class ImmunopeptidomeUIProcessing(_ImmunopeptidomeUI):
    @staticmethod
    def hla_alleles(alleles_selected: list):
        st.text(f"Selected: {sorted(alleles_selected)}")
        return alleles_selected

    @staticmethod
    def transcript_ids(transcript_ids: list):
        st.text(f"Selected: {sorted(transcript_ids)}")
        return transcript_ids

    @staticmethod
    def subset_method(transcripts: list, method: str):
        if len(transcripts) == 0 or method == "None":
            st.text("No subset selected.")
            return [], []
        elif method == "Retention":
            st.text(f"Retaining subset of transcripts: {transcripts}")
            return transcripts, []
        elif method == "Exclusion":
            st.text(f"Excluding subset of transcripts: {transcripts}")
            return [], transcripts
        else:
            raise ValueError(
                f"Method does not belong to options: "
                f"{ImmunopeptidomesUIOptions.subset_method()}"
            )

    @staticmethod
    def name(name: str):
        if not name.endswith(".bed"):
            name += ".bed"

        st.text(
            f"Output file will be saved to "
            f"{Directories.app_immunopeptidomes(name)}"
        )
        return name


class RunTab:
    @staticmethod
    def pipeline(params: Parameters):
        params.cache_dir.mkdir(exist_ok=True)
        output = st.empty()
        with st_capture(output.code):
            t_start = time()
            output = st.empty()
            with st_capture(output.code):
                run_pipeline(params)
            t_end = time()

        data_frame = pd.read_csv(params.results_path, sep="\t")

        st.text(f"Pipeline run in {int(t_end - t_start)} seconds")
        st.dataframe(data_frame, hide_index=True)
        st.text(f"dN/dS file: {params.results_path}")

    @staticmethod
    def link_vep(cache_location: pathlib.Path):
        output = st.empty()
        with st_capture(output.code):
            src_dst_links = _get_src_dst_link_pairs(cache_location)
            _link_src_dst_pairs(src_dst_links, _skip_user_input=True)

    @staticmethod
    def download(
        species: str, assembly: str, release: int, download_type: str
    ):
        if assembly == "GRCh37":
            assert species == "homo_sapiens"
            data = EnsemblData.homo_sapiens_GRCh37()
        else:
            data = EnsemblData(species=species, assembly=assembly)

        t_start = time()
        output = st.empty()
        with st_capture(output.code):
            if download_type == "toplevel":
                data.compute_all_toplevel(release)
                checks = (
                    data.toplevel_gz_done,
                    data.toplevel_fa_done,
                    data.toplevel_fai_done,
                    data.sizes_done,
                )
            else:
                data.compute_all_primary_assembly(release=release)
                checks = (
                    data.primary_assembly_gz_done,
                    data.primary_assembly_fa_done,
                    data.primary_assembly_fai_done,
                    {release},
                )

        process_ok = all([release in check for check in checks])
        st.text(f"All complete: {process_ok}")
        t_end = time()
        st.text(f"... in {int(t_end - t_start)} seconds")

    @staticmethod
    def annotate(*args, **kwargs):
        pass

    @staticmethod
    def immunopeptidome(
        hla_selections,
        output_name,
        transcripts_retained,
        transcripts_excluded,
    ):
        try:
            immunopeptidome_from_hla(
                *hla_selections,
                output_name=output_name,
                restricted_transcript_ids=transcripts_retained,
                excluded_transcript_ids=transcripts_excluded,
            )
            st.text(
                f"Completed: {Directories.app_immunopeptidomes(output_name)}"
            )
        except RuntimeError:
            st.text(
                "Process failed with currently defined options. This was "
                "likely caused by the selected HLA being unavailable in "
                "the (filtered) transcript file."
            )
