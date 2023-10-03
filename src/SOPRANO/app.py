import time

import pandas as pd
import streamlit as st

import SOPRANO.utils.path_utils
from SOPRANO.core import objects
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.app_utils import (
    get_annotated_input_options,
    get_coordinate_options,
    get_human_genome_options,
    get_immunopeptidome_options,
    st_capture,
)
from SOPRANO.utils.path_utils import Directories

genome_options = get_human_genome_options()
annotated_input_options = get_annotated_input_options()
immunopeptidome_options = get_immunopeptidome_options()
coordinate_options = get_coordinate_options()


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


if __name__ == "__main__":
    gui_tab, vep_tab, genome_tab, annotate_tab, info_tab = st.tabs(
        ["GUI", "Link VEP", "Download Genomes", "Annotator", "Information"]
    )

    with gui_tab:
        # Init app
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
