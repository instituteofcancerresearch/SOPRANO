import os
import pathlib
import shutil

import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from SOPRANO.core import objects
from SOPRANO.utils.app_utils import (
    AnnotatorUIOptions,
    AnnotatorUIProcessing,
    DownloaderUIOptions,
    DownloaderUIProcessing,
    ImmunopeptidomesUIOptions,
    ImmunopeptidomeUIProcessing,
    LinkVEPUIProcessing,
    PipelineUIOptions,
    PipelineUIProcessing,
    RunTab,
    process_text_and_file_inputs,
)
from SOPRANO.utils.path_utils import Directories


def with_tab_pipeline(tab: DeltaGenerator):
    with tab:
        st.title("Pipeline deployment")
        st.markdown(
            "To run SOPRRANO, you must choose three core inputs:\n"
            "- reference genome\n"
            "- annotated mutation file\n"
            "- bed protein file\n\n"
            "followed by the run time parameters.\n\n"
            "Once all required fields have been filled, you can deploy with"
            "'Run Pipeline'."
        )

        genome_selection = st.selectbox(
            "Select a reference genome:",
            PipelineUIOptions.genome_reference(),
            help="Can't find what you're looking for? Go to 'Step 1'.",
        )
        genome_ready, genome_processed = PipelineUIProcessing.genome_reference(
            genome_selection
        )

        annotation_selection = st.selectbox(
            "Select a VEP annotated file:",
            PipelineUIOptions.annotated_mutations(),
            help="Build your own custom files with 'Step 2'.",
        )
        (
            annotation_ready,
            annotation_processed,
        ) = PipelineUIProcessing.annotated_mutations(annotation_selection)

        immunopeptidome_selection = st.selectbox(
            "Select a BED protein file:",
            PipelineUIOptions.immunopeptidome(),
            help="Build your own custom files with 'Step 3'.",
        )
        (
            immunopeptidome_ready,
            immunopeptidome_processed,
        ) = PipelineUIProcessing.immunopeptidome(immunopeptidome_selection)

        substitution_selection = st.selectbox(
            "Trinucleotide context correction:",
            PipelineUIOptions.substitution_method(),
        )
        (
            substitution_ready,
            substitution_processed,
        ) = PipelineUIProcessing.substitution_method(substitution_selection)

        exclude_drivers = st.checkbox("Exclude driver genes:", value=True)
        use_randomization = st.checkbox("Use randomization:", value=False)

        coordinates_ready = True
        random_seed = -1
        coordinates_processed = None

        if use_randomization:
            random_seed = st.number_input(
                "Select random seed for randomization:",
                min_value=-1,
                value="min",
                help="Choose a seed value to replicate dN/dS derived from a "
                "randomized background. If '-1' is selected, then no seed "
                "value is selected, and results will not be reproducible.",
            )
            coordinates_selection = st.selectbox(
                "Coordinate file:",
                PipelineUIOptions.coordinates(),
                help="A BED coordinate file can be used to constrain the "
                "regions of randomization.",
            )
            (
                coordinates_ready,
                coordinates_processed,
            ) = PipelineUIProcessing.coordinates(coordinates_selection)

        cache_selected = st.text_input(
            "Cache directory:",
            value=Directories.soprano_cache().as_posix(),
            help="This is the root folder for caching results (it must exist)."
            " A sub-folder will be created using the name of the "
            "analysis as given below.",
        )
        cache_ready, cache_processed = PipelineUIProcessing.cache(
            cache_selected
        )

        job_name_selection = st.text_input(
            "Define a name for the output of your analysis:"
        )
        job_name_ready, job_name_processed = PipelineUIProcessing.job_name(
            job_name_selection
        )

        ready = (
            genome_ready
            and annotation_ready
            and annotation_ready
            and immunopeptidome_ready
            and substitution_ready
            and coordinates_ready
            and cache_ready
            and job_name_ready
        )

        if st.button("Run Pipeline", disabled=not ready):
            params = objects.Parameters(
                analysis_name=job_name_selection,
                input_path=annotation_processed,
                bed_path=immunopeptidome_processed,
                cache_dir=job_name_processed,
                random_regions=coordinates_processed,
                use_ssb192=substitution_processed == 192,
                use_random=use_randomization,
                exclude_drivers=exclude_drivers,
                seed=random_seed,
                transcripts=objects.TranscriptPaths.defaults(),
                genomes=genome_processed,
            )
            RunTab.pipeline(params=params)


def with_tab_genomes(tab: DeltaGenerator):
    with tab:
        st.subheader("Link existing reference genome files")
        st.markdown(
            "SOPRANO uses a self-contained directory structure for the "
            "the download and caching of derived genomic data.\n\n"
            "Users who have an existing Ensembl VEP configuration on their "
            "computer may have readily available reference genomes "
            "downloaded. These downloads can be linked to the SOPRANO "
            "directories by running the button below.\n\nNote: "
            " the default VEP cache that is searched for is `~/.vep` but you "
            "can define any other non-standard locations that reference "
            "genomes may be found within."
        )

        cache_location_selection = st.text_input(
            "VEP cache location:", value=Directories.std_sys_vep().as_posix()
        )

        cache_location_processed = LinkVEPUIProcessing.cache_location(
            cache_location_selection
        )

        if st.button("Attempt VEP cache link", disabled=False):
            RunTab.link_vep(cache_location_processed)

        st.subheader("Download new reference genome files")
        st.markdown(
            "You can use SOPRANO to download reference genomes from the "
            "ensembl FTP server by making use of the below definitions, before"
            " clicking `Download`.\n\n"
            "The core SOPRANO calculation requires toplevel reference data "
            "to be available. Secondary downloads of the primary assembly "
            "may be used to accelerate the annotation procedure; though this "
            "is _not_ essential."
        )
        feature_disabled = (
            os.environ.get("SOPRANO_DISABLE_DOWNLOADS", "False") == "True"
        )

        if feature_disabled:
            st.warning(
                "Downloading genome references is currently disabled at this "
                "app hosting.\n\n"
                "Contact the administrators scientificcomputingteam@icr.ac.uk "
                "to request additional genomic files to be added."
            )

        species_selection = st.text_input(
            "Define the species",
            value="Homo Sapiens",
            disabled=True,
        )
        species_processed = DownloaderUIProcessing.species(species_selection)

        assembly_selection = st.text_input(
            "Define the genome reference:",
            value="GRCh38",
        )
        assembly_processed = DownloaderUIProcessing.assembly(
            assembly_selection
        )

        release_selection = st.number_input(
            "Define the Ensembl release:",
            min_value=76,
            key="download_release",
            value=110,
        )
        release_processed = DownloaderUIProcessing.release(release_selection)

        type_selection = st.selectbox(
            "Download type:",
            options=DownloaderUIOptions.type(),
        )
        type_processed = DownloaderUIProcessing.type(type_selection)

        if st.button("Download", disabled=feature_disabled):
            RunTab.download(
                species=species_processed,
                assembly=assembly_processed,
                release=release_processed,
                download_type=type_processed,
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


def with_tab_annotator(tab: DeltaGenerator):
    with tab:
        anno_cache = AnnoCache()
        tmp_dir = anno_cache.path

        st.title("Annotate VCF files")

        st.markdown("Generate an annotated mutation file from VCF files.")

        vcf_definition_method_selection = st.radio(
            "Method for defining VCF files to annoatate:",
            options=AnnotatorUIOptions.vcf_definition_method(),
        )

        if vcf_definition_method_selection == "File uploader":
            vcf_upload_selection = st.file_uploader(
                "Upload VCF files to annotate.",
                accept_multiple_files=True,
                type=["vcf", "vcf.gz"],
            )

            (
                vcf_dir_ready,
                vcf_dir_processed,
            ) = AnnotatorUIProcessing.vcf_upload_sources(
                vcf_upload_selection=vcf_upload_selection, tmp_dir=tmp_dir
            )

        else:
            vcf_dir_selection = st.text_input(
                "Directory containing VCF files:", value=pathlib.Path.home()
            )
            (
                vcf_dir_ready,
                vcf_dir_processed,
            ) = AnnotatorUIProcessing.vcf_path_sources(vcf_dir_selection)

        assembly_selection = st.selectbox(
            "Genome reference assembly:",
            options=AnnotatorUIOptions.assembly_type(),
        )
        (
            assembly_ready,
            assembly_processed,
        ) = AnnotatorUIProcessing.assembly_type(assembly_selection)

        st.markdown(
            "Note: _For multiple file uploads, this name will be applied "
            "to the concatenation of output annotated files. "
            "Otherwise no concatenation will be generated._"
        )

        name_selection = st.text_input(
            "Choose a name for the annotated output:"
        )

        name_ready, name_processed = AnnotatorUIProcessing.output_name(
            name_selection
        )

        if st.button(
            "Annotate",
            disabled=not (vcf_dir_ready and assembly_ready and name_ready),
        ):
            RunTab.annotate(
                source_path=vcf_dir_processed,
                output_name=name_processed,
                assembly=assembly_processed,
            )
            st.text(f"Processed sources @ {vcf_dir_processed}")
            anno_cache.clean_up()


def with_tab_info(tab: DeltaGenerator):
    with tab:
        st.title("Welcome to SOPRANO! :wave:")
        st.caption("Selection On PRotein ANnotated regiOns")
        st.markdown(
            "This application is designed to provide a user interface to the "
            "SOPRANO computational pipeline, without the need of command line "
            "intervention."
            "\n\n"
            "There are three essential files required to run "
            "SOPRANO. These define the\n"
            "1. Reference genome\n"
            "2. Annotated somatic mutations\n"
            "3. Immunopeptidome\n"
            "\n\n"
            "These three inputs can be configured in term via the tabs "
            "indicating steps 1, 2 and 3. Once you have prepared your data, "
            "step 4 will enable you to run the pipeline, subject to "
            "further runtime configuration choices."
            "\n\n"
            "Any technical issues can be raised on [GitHub]"
            "(https://github.com/instituteofcancerresearch/SOPRANO/issues)"
        )


def with_tab_immunopeptidome(tab: DeltaGenerator):
    with tab:
        st.header("HLA-allele selections")

        st.markdown(
            "In order to generate custom, or, patient specific "
            "immunopeptidome files to run through SOPRANO, we can "
            "specify a set of HLA alleles to restrict the analysis to."
        )

        hla_alleles_selected = st.multiselect(
            "Select 1-6 HLA alleles:",
            options=ImmunopeptidomesUIOptions.hla_alleles(),
        )

        (
            hla_alleles_ready,
            hla_alleles_processed,
        ) = ImmunopeptidomeUIProcessing.hla_alleles(hla_alleles_selected)

        filter_by_transcript = st.checkbox(
            "Filter by ensembl transcript IDs?", value=False
        )

        if filter_by_transcript:
            st.header("Ensembl transcript selections")

            transcript_supply_method = st.radio(
                "Select a method to define transcript IDs:",
                options=("File uploader", "Text box"),
            )

            if transcript_supply_method == "File uploader":
                raw_transcript_input = st.file_uploader(
                    "Upload a file containing ensembl "
                    "transcript IDs over separate lines: ",
                )
            else:
                _raw_transcript_input = st.text_area(
                    "Enter ensembl transcript IDs " "over separate lines:",
                    value="",
                )

                if _raw_transcript_input:
                    raw_transcript_input = _raw_transcript_input
                else:
                    raw_transcript_input = None

            (
                processed_transcripts_ready,
                processed_transcript_input,
            ) = process_text_and_file_inputs(raw_transcript_input)

            subset_method_selected = st.radio(
                "Select method to subset available Ensembl transcript IDs:",
                options=ImmunopeptidomesUIOptions.subset_method(),
            )

            (
                subset_ready,
                retained_excluded,
            ) = ImmunopeptidomeUIProcessing.subset_method(
                processed_transcript_input, subset_method_selected
            )

        else:
            processed_transcripts_ready = True
            subset_ready = True
            retained_excluded = [], []

        transcripts_ready = processed_transcripts_ready and subset_ready
        transcripts_retained, transcripts_excluded = retained_excluded

        name_selected = st.text_input("Immunopeptidome name:")
        name_ready, name_processed = ImmunopeptidomeUIProcessing.name(
            name_selected
        )

        selections_ready = hla_alleles_ready and transcripts_ready

        if not selections_ready:
            st.warning(
                "Please provide a set of alleles (and transcript IDs if "
                "relevant) to construct the restricted immunopeptidome "
                "input file."
            )

        if not name_ready:
            st.warning("Please provide an input file name.")

        if st.button(
            "Build immunopeptidome",
            disabled=not (selections_ready and name_ready),
        ):
            RunTab.immunopeptidome(
                hla_alleles_processed,
                output_name=name_processed,
                transcripts_retained=transcripts_retained,
                transcripts_excluded=transcripts_excluded,
            )


if __name__ == "__main__":
    st.set_page_config(layout="wide")
    (
        welcome_tab,
        genome_tab,
        annotate_tab,
        immunopeptidome_tab,
        pipeline_tab,
    ) = st.tabs(
        [
            "Welcome!",
            "Step 1: Prepare genome reference",
            "Step 2: Annotate mutations",
            "Step 3: Prepare immunopeptidome",
            "Step 4: Run pipeline",
        ]
    )
    with_tab_info(welcome_tab)
    with_tab_genomes(genome_tab)
    with_tab_annotator(annotate_tab)
    with_tab_immunopeptidome(immunopeptidome_tab)
    with_tab_pipeline(pipeline_tab)
