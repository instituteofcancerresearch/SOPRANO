import pathlib
import time

import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

import SOPRANO.utils.path_utils
from SOPRANO.core import objects
from SOPRANO.hla2ip import (
    get_hla_options,
    get_transcript_id_options,
    immunopeptidome_from_hla,
)
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.app_utils import (
    get_annotated_input_options,
    get_coordinate_options,
    get_human_genome_options,
    get_immunopeptidome_options,
    st_capture,
)
from SOPRANO.utils.parse_utils import fix_species_arg
from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.vep_utils import (
    _get_src_dst_link_pairs,
    _link_src_dst_pairs,
)

genome_options = get_human_genome_options()
annotated_input_options = get_annotated_input_options()
immunopeptidome_options = get_immunopeptidome_options()
coordinate_options = get_coordinate_options()
hla_options = [""] + get_hla_options()
transcript_id_options = get_transcript_id_options()


def process_genome_selection():
    genome_selection = st.session_state.genome_selection

    ref_id, rel_id = genome_selection.split(" - Ensembl release ")

    (
        genomes_path,
        chroms_path,
    ) = SOPRANO.utils.path_utils.genome_pars_to_paths(ref_id, rel_id)

    st.session_state.ref_genome = objects.GenomePaths(
        sizes=chroms_path, fasta=genomes_path
    )

    st.text(f"Selected: {st.session_state.ref_genome.fasta}")


def process_annotated_input_selection():
    input_selection = st.session_state.input_selection
    st.session_state.input_path = annotated_input_options[input_selection]
    st.text(f"Selected: {st.session_state.input_path}")


def process_immunopeptidome_selection():
    bed_selection = st.session_state.bed_selection
    st.session_state.bed_path = immunopeptidome_options[bed_selection]
    st.text(f"Selected: {st.session_state.bed_path}")


def process_substitution_selection():
    subs_selection = st.session_state.subs_selection
    st.session_state.use_ssb192 = subs_selection == 192
    st.text(f"Using SSB192: {st.session_state.use_ssb192}")


def process_job_name_selection():
    job_name = st.session_state.job_name
    st.session_state.cache_dir = Directories.cache(job_name)
    st.text(f"Cache location: {st.session_state.cache_dir}")


def process_randomization_region():
    region_selection = st.session_state.random_regions
    st.session_state.random_regions_path = coordinate_options[region_selection]
    st.text(f"Selected: {st.session_state.random_regions_path}")


def process_vep_cache_selection():
    vep_cache_selection = st.session_state.vep_cache
    if isinstance(vep_cache_selection, str):
        st.session_state.vep_cache_dir = pathlib.Path(vep_cache_selection)
    else:
        st.session_state.vep_cache_dir = vep_cache_selection

    st.text(f"Selected: {st.session_state.vep_cache_dir.as_posix()}")


def process_download_species_selection():
    download_species_selection = st.session_state.download_species_selection
    st.session_state.download_species = fix_species_arg(
        download_species_selection
    )
    st.text(f"Selected: {st.session_state.download_species}")


def process_vcf_selection():
    pass  # TODO


def run_pipeline_in_app():
    st.session_state.cache_dir.mkdir(exist_ok=True)

    t_start = time.time()
    output = st.empty()
    with st_capture(output.code):
        run_pipeline(st.session_state.params)
        # run_local_ssb_selection.main(st.session_state.namespace)
    t_end = time.time()

    st.session_state.compute_time_str = (
        f"Pipeline run in {int(t_end - t_start)} seconds"
    )

    st.session_state.data_frame = pd.read_csv(
        st.session_state.params.results_path, sep="\t"
    )

    st.session_state.job_complete = True

    st.text(st.session_state.compute_time_str)
    st.dataframe(st.session_state.data_frame, hide_index=True)
    st.text(f"dN/dS file: {st.session_state.params.results_path}")


def with_tab_pipeline(tab: DeltaGenerator):
    with tab:
        st.title("SOPRANO")
        st.caption("Selection On PRotein ANnotated regiOns")

        # Derived genome definitions
        st.selectbox(
            "Select a reference genome:",
            genome_options.keys(),
            key="genome_selection",
        )
        process_genome_selection()

        # VEP annotated file
        st.selectbox(
            "Select a VEP annotated file:",
            annotated_input_options.keys(),
            key="input_selection",
        )
        process_annotated_input_selection()

        # BED protein transcript file
        st.selectbox(
            "Select a BED protein file:",
            immunopeptidome_options.keys(),
            key="bed_selection",
        )
        process_immunopeptidome_selection()

        # Substitution method
        st.selectbox(
            "Select a substitution method:", (192, 7), key="subs_selection"
        )
        process_substitution_selection()

        # Exclude driver genes
        st.checkbox("Exclude driver genes:", value=True, key="exclude_drivers")

        # Use randomization
        st.checkbox("Use randomization:", value=False, key="use_random")

        if st.session_state.use_random:
            # Select random seed
            st.number_input(
                "Select random seed for randomization:",
                min_value=-1,
                value="min",
                key="random_seed",
            )
            st.selectbox(
                "Select a BED file defining the regions to randomize over:",
                coordinate_options.keys(),
                key="random_regions",
            )
            process_randomization_region()
        else:
            st.session_state.random_seed = -1
            st.session_state.random_regions_path = None

        # Pipeline job name & cache
        st.text_input(
            "Define a name for the output of your analysis:", key="job_name"
        )
        process_job_name_selection()

        st.session_state.params = objects.Parameters(
            analysis_name=st.session_state.job_name,
            input_path=st.session_state.input_path,
            bed_path=st.session_state.bed_path,
            cache_dir=st.session_state.cache_dir,
            random_regions=st.session_state.random_regions_path,
            use_ssb192=st.session_state.use_ssb192,
            use_random=st.session_state.use_random,
            exclude_drivers=st.session_state.exclude_drivers,
            seed=st.session_state.random_seed,
            transcripts=objects.TranscriptPaths.defaults(),
            genomes=st.session_state.ref_genome,
        )

        st.session_state.job_complete = False

        if st.button("Run Pipeline"):
            run_pipeline_in_app()


def with_tab_vep(tab: DeltaGenerator):
    with tab:
        st.title("Link VEP")
        st.caption(
            f"Symbolically link files in your VEP cache to the SOPRANO data "
            f"folder: {Directories.data()}"
        )

        default_location = Directories.std_sys_vep()

        st.text_input(
            "VEP cache location:", key="vep_cache", value=default_location
        )
        process_vep_cache_selection()

        if st.button("Link", disabled=False):
            output = st.empty()
            with st_capture(output.code):
                src_dst_links = _get_src_dst_link_pairs(
                    st.session_state.vep_cache_dir
                )
                _link_src_dst_pairs(src_dst_links, _skip_user_input=True)


def with_tab_genomes(tab: DeltaGenerator):
    with tab:
        st.title("Download Reference Genomes")
        st.caption(
            "Download Ensembl genome reference data. See "
            "https://www.ensembl.org/info/genome/variation/species/"
            "species_data_types.html for definitions."
        )

        st.text_input(
            "Define the species",
            value="Homo Sapiens",
            key="download_species_selection",
        )
        process_download_species_selection()

        st.text_input(
            "Define the genome reference:",
            value="GRCh38",
            key="download_assembly",
        )
        st.text(f"Selected: {st.session_state.download_assembly}")

        st.number_input(
            "Define the Ensembl release:",
            min_value=76,
            key="download_release",
            value=110,
        )

        msg = f"Selected: {st.session_state.download_release}"
        if st.session_state.download_release > 110:
            msg += "\n\n[Warning Oct 1 2023: Latest Ensembl release is 110]"
        st.text(msg)

        st.selectbox(
            "Download type:",
            options=("toplevel", "primary_assembly"),
            key="download_type",
        )
        st.text(f"Selected: {st.session_state.download_type}")

        if st.button("Download", disabled=False):
            if st.session_state.download_assembly == "GRCh37":
                assert st.session_state.download_species == "homo_sapiens"
                data = objects.EnsemblData.homo_sapiens_GRCh37()
            else:
                data = objects.EnsemblData(
                    species=st.session_state.download_species,
                    assembly=st.session_state.download_assembly,
                )

            t_start = time.time()
            output = st.empty()
            with st_capture(output.code):
                if st.session_state.download_type == "toplevel":
                    data.compute_all_toplevel(
                        st.session_state.download_release
                    )
                    checks = (
                        data.toplevel_gz_done,
                        data.toplevel_fa_done,
                        data.toplevel_fai_done,
                        data.sizes_done,
                    )
                else:
                    data.compute_all_primary_assembly(
                        st.session_state.download_release
                    )
                    checks = (
                        data.primary_assembly_gz_done,
                        data.primary_assembly_fa_done,
                        data.primary_assembly_fai_done,
                        {st.session_state.download_release},
                    )
                process_ok = all(
                    [
                        st.session_state.download_release in check
                        for check in checks
                    ]
                )
                print(f"All downloads complete: {process_ok}")
            t_end = time.time()
            st.text(f"Download complete in {int(t_end - t_start)} seconds")


def with_tab_annotator(tab: DeltaGenerator):
    with tab:
        st.title("Annotate VCF File")
        st.caption(
            "Upload a VCF file to annotate for use in the SOPRANO pipeline."
        )
        st.file_uploader("Select a VCF file:", key="vcf_selection")
        process_vcf_selection()
        if st.button("Annotate", disabled=True):
            pass  # TODO


def with_tab_info(tab: DeltaGenerator):
    with tab:
        st.title("Information")
        st.caption("Description of what is going on...")


def with_tab_immunopeptidome(tab: DeltaGenerator):
    with tab:
        st.title("Immunopeptidomes")
        st.caption(
            "Generate a patient specific immunopeptidome file derived from "
            "the SOPRANO internal master file and HLA allele selections."
        )

        st.multiselect(
            "Select HLA alleles:",
            options=hla_options,
            key="hla_selections",
            max_selections=4,
        )

        st.text_input("Immunopeptidome name:", key="input_ip_name")

        if not st.session_state.input_ip_name.endswith(".bed"):
            st.session_state.ip_name = st.session_state.input_ip_name + ".bed"
        else:
            st.session_state.ip_name = st.session_state.input_ip_name

        st.text(
            f"Output path will be: "
            f"{Directories.app_immunopeptidomes(st.session_state.ip_name)}"
        )

        _no = "No"
        _retain = "Retain subset of IDs"
        _exclude = "Exclude subset of IDs"

        st.selectbox(
            "Filter transcript IDs:",
            options=(_no, _retain, _exclude),
            key="filter_by_transcript",
        )

        if st.session_state.filter_by_transcript != _no:
            st.multiselect(
                "Select transcript IDs:",
                options=transcript_id_options,
                key="selected_transcripts",
            )

            if st.session_state.filter_by_transcript == _retain:
                st.session_state.transcripts_retained = (
                    st.session_state.selected_transcripts
                )
                st.session_state.transcripts_excluded = []
            else:
                st.session_state.transcripts_retained = []
                st.session_state.transcripts_excluded = (
                    st.session_state.selected_transcripts
                )
        else:
            st.session_state.selected_transcripts = []
            st.session_state.transcripts_retained = []
            st.session_state.transcripts_excluded = []

        # st.text(st.session_state.selected_transcripts)
        # st.text(st.session_state.transcripts_retained)
        # st.text(st.session_state.transcripts_excluded)

        if st.button("Build immunopeptidome"):
            immunopeptidome_from_hla(
                *st.session_state.hla_selections,
                output_name=st.session_state.ip_name,
                restricted_transcript_ids=st.session_state.transcripts_retained,
                excluded_transcript_ids=st.session_state.transcripts_excluded,
            )
            st.text(f"Completed: {st.session_state.ip_name}")


if __name__ == "__main__":
    st.set_page_config(layout="wide")
    (
        pipeline_tab,
        vep_tab,
        genome_tab,
        annotate_tab,
        immunopeptidome_tab,
        info_tab,
    ) = st.tabs(
        [
            "Pipeline",
            "Link VEP",
            "Download Genomes",
            "Annotator",
            "Immunopeptidomes",
            "Information",
        ]
    )

    with_tab_pipeline(pipeline_tab)
    with_tab_vep(vep_tab)
    with_tab_genomes(genome_tab)
    with_tab_annotator(annotate_tab)
    with_tab_immunopeptidome(immunopeptidome_tab)
    with_tab_info(info_tab)
