import os
import pathlib
import shutil
from contextlib import contextmanager, redirect_stdout
from io import StringIO
from time import time

import pandas as pd
import streamlit as st
from streamlit.runtime.uploaded_file_manager import UploadedFile

from SOPRANO.core.objects import EnsemblData, Parameters
from SOPRANO.hla2ip import immunopeptidome_from_hla
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils import anno_utils
from SOPRANO.utils.parse_utils import fix_species_arg
from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe
from SOPRANO.utils.vep_utils import (
    _get_src_dst_link_pairs,
    _link_src_dst_pairs,
)


def _lines_ok(lines: list | tuple, min_args: int, max_args: int):
    if min_args > max_args:
        raise ValueError(f"(min = {min_args}) > (max = {max_args})")

    return min_args <= len(lines) <= max_args


def process_text_and_file_inputs(
    raw_input: str | UploadedFile | None,
    min_args=0,
    max_args=int(1e9),
    remove_empty_lines=True,
):
    if raw_input is None:
        return False, []
    elif isinstance(raw_input, str):
        if raw_input == "":
            return False, []
        else:
            lines = raw_input.split("\n")
    elif isinstance(raw_input, UploadedFile):
        lines = StringIO(raw_input.getvalue().decode("utf-8")).readlines()
    else:
        raise TypeError(raw_input)

    lines = [line.strip() for line in lines]

    if remove_empty_lines:
        lines = [line for line in lines if line != ""]

    status_ok = _lines_ok(lines, min_args, max_args)

    if not status_ok:
        st.warning(
            f"Number of arguments parsed from input not within bounds "
            f"[{min_args}, {max_args}]"
        )

    return status_ok, lines


def text_or_file(
    desc: str,
    min_args: int = 0,
    max_args: int = int(1e9),
    help_text: str | None = None,
    help_upload: str | None = None,
):
    raw_text_input = st.text_area(desc, value="", help=help_text)
    raw_file_input = st.file_uploader(desc, help=help_upload)

    text_ready, text_input = process_text_and_file_inputs(
        raw_text_input, min_args=min_args, max_args=max_args
    )
    file_ready, file_input = process_text_and_file_inputs(
        raw_file_input, min_args=min_args, max_args=max_args
    )

    if text_ready == file_ready:
        ready = False
        content = None

        if text_ready:
            st.warning(
                "Multiple input selections detected!"
                " Provide manual text input OR upload a file."
            )

    elif text_ready:
        ready = True
        content = text_input
    elif file_ready:
        ready = True
        content = file_input
    else:
        ready = False
        content = None

    assert isinstance(content, list | None), (content, type(content))

    # Ready status should disable button prompt in UI
    return ready, content


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

    @staticmethod
    def job_name(*args, **kwargs):
        pass

    @staticmethod
    def cache(*args, **kwargs):
        pass


class PipelineUIOptions(_PipelineUI):
    OPTION_SUBS_METHOD_SSB7 = "SSB7"
    OPTION_SUBS_METHOD_SSB192 = "SSB192"

    @staticmethod
    def genome_reference():
        homo_sapiens_dir = Directories.homo_sapien_reference_files()

        genome_dirs = [
            item for item in homo_sapiens_dir.glob("*") if item.is_dir()
        ][::-1]

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

        # Note: Dictionary keys are presented as options in the UI!

        return options_dict

    @staticmethod
    def annotated_mutations():
        options_dict = {}
        for directory in (
            Directories.annotation_example_files(),
            Directories.app_annotated_inputs(),
        ):
            for x in directory.glob("*.anno*"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def immunopeptidome():
        options_dict = {}
        for directory in (
            Directories.immunopeptidome_example_files(),
            Directories.app_immunopeptidomes(),
        ):
            for x in directory.glob("*.bed"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def substitution_method():
        return (
            PipelineUIOptions.OPTION_SUBS_METHOD_SSB192,
            PipelineUIOptions.OPTION_SUBS_METHOD_SSB7,
        )

    @staticmethod
    def coordinates():
        options_dict = {}
        for x in Directories.app_coordinate_files().glob("*.bed"):
            options_dict[x.name] = x
        return options_dict


class PipelineUIProcessing(_PipelineUI):
    @staticmethod
    def genome_reference(genome_selection: str | None):
        if genome_selection is None:
            st.warning("Warning: No genome selection.")
            genome_ready = False
            genome_reference_path = None
        else:
            genome_ready = True
            assembly, release = genome_selection.split(" - Ensembl release ")
            data = EnsemblData(species="homo_sapiens", assembly=assembly)
            fasta_path = data.toplevel_fa_path(int(release))
            chrom_path = data.toplevel_chrom_path(int(release))
            st.text(
                f"Selected reference: {fasta_path}\n"
                f"Selected chrom sizes: {chrom_path}"
            )
            genome_reference_path = data.get_genome_reference_paths(
                int(release)
            )

        return genome_ready, genome_reference_path

    @staticmethod
    def annotated_mutations(annotation_selection: str):
        options_dict = PipelineUIOptions.annotated_mutations()
        return True, options_dict[annotation_selection]

    @staticmethod
    def immunopeptidome(immunopeptidome_selection: str):
        options_dict = PipelineUIOptions.immunopeptidome()
        return True, options_dict[immunopeptidome_selection]

    @staticmethod
    def substitution_method(subs_selection: str):
        if subs_selection == PipelineUIOptions.OPTION_SUBS_METHOD_SSB192:
            return True, 192
        elif subs_selection == PipelineUIOptions.OPTION_SUBS_METHOD_SSB7:
            return True, 7
        else:
            raise ValueError(f"Unrecognized method: {subs_selection}")

    @staticmethod
    def coordinates(coordinates_selection: str):
        options_dict = PipelineUIOptions.coordinates()
        return True, options_dict.get(coordinates_selection, None)

    @staticmethod
    def job_name(job_name: str):
        if job_name == "":
            job_name_ready = False
            cache_dir = None
        else:
            cache_dir = Directories.soprano_cache(job_name)

            if cache_dir.exists():
                st.warning(f"Pipeline cache already in use: {cache_dir}")
                job_name_ready = False
            else:
                st.text(f"Pipeline results cache: {cache_dir}")
                job_name_ready = True

        return job_name_ready, cache_dir

    @staticmethod
    def cache(cache_selected: str):
        if os.path.exists(cache_selected):
            os.environ["SOPRANO_CACHE"] = cache_selected
            # st.text(f"Selected: {cache_selected}")
            cache_ready = True
        else:
            st.warning(f"Cache directory does not exist: {cache_selected}")
            cache_ready = False

        return cache_ready, cache_selected


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

    @staticmethod
    def vep_cache_location(*args, **kwargs):
        pass


class DownloaderUIOptions(_DownloaderUI):
    OPTION_TYPE_TOPLEVEL = "toplevel"
    OPTION_TYPE_PRIMARY_ASSEMBLY = "primary_assembly"

    @staticmethod
    def type():
        return (
            DownloaderUIOptions.OPTION_TYPE_TOPLEVEL,
            DownloaderUIOptions.OPTION_TYPE_PRIMARY_ASSEMBLY,
        )


class DownloaderUIProcessing(_DownloaderUI):
    @staticmethod
    def species(species_selection: str):
        output = fix_species_arg(species_selection)
        st.text(f"Selected: {output}")
        return True, output

    @staticmethod
    def assembly(assembly_selection: str):
        st.text(f"Selected: {assembly_selection}")
        return True, assembly_selection

    @staticmethod
    def release(release: str):
        output = int(release)
        st.text(f"Selected: {output}")

        if output > 110:
            st.text("[Warning] Oct 1 2023: Latest Ensembl release is 110")

        return True, output

    @staticmethod
    def type(type_selection: str):
        if type_selection not in (
            DownloaderUIOptions.OPTION_TYPE_TOPLEVEL,
            DownloaderUIOptions.OPTION_TYPE_PRIMARY_ASSEMBLY,
        ):
            raise ValueError(type_selection)

        st.text(f"Selected: {type_selection}")
        return True, type_selection

    @staticmethod
    def vep_cache_location(cache_location: str):
        output = pathlib.Path(cache_location)

        ready = output.exists()

        if not ready:
            st.warning(f"Cache directory does not exist: {cache_location}")
        else:
            st.text(f"Selected: {output}")
        return ready, output


class _AnnotatorUI:
    @staticmethod
    def vcf_definition_method(*args, **kwargs):
        pass

    @staticmethod
    def vcf_upload_sources(*args, **kwargs):
        pass

    @staticmethod
    def vcf_path_sources(*args, **kwargs):
        pass

    @staticmethod
    def assembly_type(*args, **kwargs):
        pass

    @staticmethod
    def output_name(*args, **kwargs):
        pass


class AnnotatorUIOptions(_AnnotatorUI):
    OPTION_VCF_DEF_METHOD_SYS_PATH = "System path"
    OPTION_VCF_DEF_METHOD_UPLOADER = "File uploader"

    OPTION_ASSEMBLY_TYPE_38 = "GRCh38"
    OPTION_ASSEMBLY_TYPE_37 = "GRCh37"

    @staticmethod
    def vcf_definition_method():
        return (
            AnnotatorUIOptions.OPTION_VCF_DEF_METHOD_SYS_PATH,
            AnnotatorUIOptions.OPTION_VCF_DEF_METHOD_UPLOADER,
        )

    @staticmethod
    def assembly_type():
        return (
            AnnotatorUIOptions.OPTION_ASSEMBLY_TYPE_38,
            AnnotatorUIOptions.OPTION_ASSEMBLY_TYPE_37,
        )


class AnnotatorUIProcessing(_AnnotatorUI):
    @staticmethod
    def vcf_upload_sources(
        vcf_upload_selection: list[UploadedFile], tmp_dir: pathlib.Path
    ):
        if not tmp_dir.exists():
            raise NotADirectoryError(tmp_dir)

        n_uploads = len(vcf_upload_selection)

        if n_uploads < 1:
            upload_ready = False
            st.warning("No vcf files uploaded to annotate.")
        else:
            upload_ready = True
            st.text(f"{n_uploads} vcf files will be processed.")

        for file in vcf_upload_selection:
            bytes = file.getvalue()
            filename = file.name
            cached_path = tmp_dir / filename

            with open(cached_path, "wb") as bytes_file:
                bytes_file.write(bytes)

        return upload_ready, tmp_dir

    @staticmethod
    def assembly_type(genome_assembly_selection: str):
        if genome_assembly_selection in (
            AnnotatorUIOptions.OPTION_ASSEMBLY_TYPE_38,
            AnnotatorUIOptions.OPTION_ASSEMBLY_TYPE_37,
        ):
            assembly_ready = True
        else:
            st.warning("Currently only supporting GRCh38 and GRCh37.")
            assembly_ready = False

        return assembly_ready, genome_assembly_selection

    @staticmethod
    def vcf_path_sources(sources_path_selection: str):
        sources_path = pathlib.Path(sources_path_selection)

        try:
            detected = anno_utils.find_vcf_files(sources_path)
        except anno_utils.NoVCFs:
            detected = []

        detected = [d.as_posix() for d in detected]

        if len(detected) < 1:
            vcf_files_ready = False
            st.warning(f"No vcf files detected in {sources_path}")
        else:
            vcf_files_ready = True
            vcf_text = "\n".join(detected)
            st.text(f"Annotated file will be constructed from: \n{vcf_text}")

        return vcf_files_ready, sources_path

    @staticmethod
    def output_name(name: str):
        if name != "":
            dest = Directories.app_annotated_inputs(name).with_suffix(
                ".vcf.anno"
            )
            st.text(f"Output path: {dest}")
        else:
            st.warning(
                "Name will be automatically constructed using input files."
            )

        name_ready = True
        return name_ready, name


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
    OPTION_SUBSET_METHOD_RETAIN = "Retain only the chosen transcript IDs"
    OPTION_SUBSET_METHOD_EXCLUDE = "Exclude all the chosen transcript IDs"

    @staticmethod
    def hla_alleles():
        hla_types_path = Directories.immunopeptidome_aux_files(
            "TCGA_hlaTypesAll.tsv"
        )
        options = pipe(
            ["cut", "-f3", hla_types_path.as_posix()],
            ["tr", ",", "\n"],
            ["sort", "-u"],
        ).split("\n")
        return options

    @staticmethod
    def transcript_ids():
        hla_binders_path = Directories.immunopeptidome_aux_files(
            "allhlaBinders_exprmean1.IEDBpeps.bed.unique_ids"
        )

        with open(hla_binders_path, "r") as f:
            transcript_options = f.read()
        # eol marker generates empty so excluded
        return transcript_options.split("\n")[:-1]

    @staticmethod
    def subset_method():
        return (
            ImmunopeptidomesUIOptions.OPTION_SUBSET_METHOD_RETAIN,
            ImmunopeptidomesUIOptions.OPTION_SUBSET_METHOD_EXCLUDE,
        )


class ImmunopeptidomeUIProcessing(_ImmunopeptidomeUI):
    @staticmethod
    def hla_alleles(alleles_selected: list):
        st.text(f"Selected: {sorted(alleles_selected)}")
        alleles_ready = len(alleles_selected) > 0
        return alleles_ready, alleles_selected

    @staticmethod
    def transcript_ids(transcript_ids: list | None):
        if transcript_ids is None:
            st.text("No transcript IDs selected.")
            return True, []
        else:
            st.text(f"Selected: {sorted(transcript_ids)}")
            ready = len(transcript_ids) > 0
            return ready, transcript_ids

    @staticmethod
    def subset_method(transcripts: list[str], method: str):
        retained: list[str] = []
        excluded: list[str] = []
        if len(transcripts) < 1:
            st.warning("No transcripts currently defined.")
            ready = False
        elif method == ImmunopeptidomesUIOptions.OPTION_SUBSET_METHOD_RETAIN:
            st.text(f"Retaining subset of transcripts: {transcripts}")
            ready = True
            retained = transcripts
        elif method == ImmunopeptidomesUIOptions.OPTION_SUBSET_METHOD_EXCLUDE:
            st.text(f"Excluding subset of transcripts: {transcripts}")
            ready = True
            excluded = transcripts
        else:
            raise ValueError(
                f"Method does not belong to options: "
                f"{ImmunopeptidomesUIOptions.subset_method()}"
            )

        retained_excluded = retained, excluded

        return ready, retained_excluded

    @staticmethod
    def name(name: str):
        if not name.endswith(".bed"):
            name += ".bed"

        output_path = Directories.app_immunopeptidomes(name)

        st.text(f"Output file will be saved to {output_path}")

        if name == ".bed":
            st.warning("Immunopeptidome output name is required.")
            ready = False
        elif output_path.exists():
            st.warning(
                f"Output file already exists: {output_path}\n"
                f"Please change or delete existing file"
            )
            ready = False
        else:
            ready = True

        return ready, name


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
    def annotate(
        source_path: pathlib.Path,
        output_name: str,
        assembly: str,
    ):
        running_msg = st.warning(
            "Annotation in progress ... please wait until this "
            "process has finished."
        )
        all_annotated_paths = anno_utils.annotate_source(
            source_path=source_path,
            output_name=None if output_name == "" else output_name,
            cache_directory=Directories.app_annotated_inputs(),
            assembly=assembly,
        )

        running_msg.empty()

        n_header_lines = 5

        for path in all_annotated_paths:
            if path.exists():
                st.success(f"Successful annotation: {path}")

                head_lines = pd.read_csv(
                    path, delimiter="\t", header=None
                ).head(n_header_lines)

                st.dataframe(head_lines)

                with open(path, "r") as f:
                    st.download_button(
                        "Download Full Annotation", f, file_name=path.name
                    )

            else:
                st.error(f"Failed annotation: {path}")

    @staticmethod
    def immunopeptidome(
        hla_selections,
        output_name,
        transcripts_retained,
        transcripts_excluded,
    ):
        n_header_lines = 5

        running_msg = st.warning(
            "Immunopeptidome processing in progress ... "
            "please wait until this process has finished."
        )
        try:
            path = immunopeptidome_from_hla(
                *hla_selections,
                output_name=output_name,
                restricted_transcript_ids=transcripts_retained,
                excluded_transcript_ids=transcripts_excluded,
            )
            st.text(
                f"Completed: {Directories.app_immunopeptidomes(output_name)}"
            )
            running_msg.empty()
            st.success(f"Processed successfully: {path}")

            head_lines = pd.read_csv(path, delimiter="\t", header=None).head(
                n_header_lines
            )

            st.dataframe(head_lines)

            with open(path, "r") as f:
                st.download_button(
                    "Download Full Immunopeptidome", f, file_name=path.name
                )

        except RuntimeError:
            st.error(
                "Process failed with currently defined options. This was "
                "likely caused by the selected HLA being unavailable in "
                "the (filtered) transcript file."
            )


class AnnoCache:
    def __init__(self):
        self.path = pathlib.Path("/tmp") / "soprano-anno"
        try:
            self.clean_up()
        except FileNotFoundError:
            pass
        finally:
            self.path.mkdir(exist_ok=True)

    def clean_up(self):
        shutil.rmtree(self.path)
