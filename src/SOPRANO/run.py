import subprocess

from SOPRANO import hla2ip
from SOPRANO.core import objects
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.anno_utils import annotate_source
from SOPRANO.utils.mpi_utils import RANK
from SOPRANO.utils.parse_utils import (
    parse_args,
    parse_genome_args,
    parse_hla,
    parse_link_vep_cache_args,
    parse_vcf_sources,
)
from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.print_utils import startup_output
from SOPRANO.utils.vep_utils import (
    _get_src_dst_link_pairs,
    _link_src_dst_pairs,
)


def run_cli(_namespace=None):
    cli_args = parse_args() if _namespace is None else _namespace
    global_parameters = objects.GlobalParameters.from_namespace(cli_args)
    startup_output(**global_parameters.__dict__)

    worker_samples_parameters = global_parameters.get_worker_samples()

    for parameters in worker_samples_parameters:
        print(
            f"Worker %05d running pipeline component: "
            f"{parameters.analysis_name}" % RANK
        )

        run_pipeline(parameters)

    global_parameters.gather()


def run_app():
    app_path = Directories.src_root("app.py")
    subprocess.run(["streamlit", "run", app_path.as_posix()])


def link_vep_cache():
    src_cache = parse_link_vep_cache_args()
    src_dst_links = _get_src_dst_link_pairs(src_cache)
    _link_src_dst_pairs(src_dst_links)
    # TODO: Just use _link_vep_cache !


def download_genome():
    args = parse_genome_args()
    startup_output(**args.__dict__)
    print("Starting download...")

    if args.assembly == "GRCh37":
        assert args.species == "homo_sapiens"
        ensembl_data = objects.EnsemblData.homo_sapiens_GRCh37()
    else:
        ensembl_data = objects.EnsemblData(args.species, args.assembly)

    if args.primary_assembly:
        ensembl_data.download_primary_assembly(args.release)
        if not args.download_only:
            ensembl_data.compute_all_primary_assembly(args.release)
    else:
        ensembl_data.download_toplevel(args.release)
        if not args.download_only:
            ensembl_data.compute_all_toplevel(args.release)


def hla2pip():
    args = parse_hla()
    hla2ip.immunopeptidome_from_hla(
        *args.hla_values,
        output_name=args.output_id,
        cache_loc=args.cache_dir,
        restricted_transcript_ids=args.restricted_transcript_ids,
        excluded_transcript_ids=args.excluded_transcript_ids,
    )


def annotate_vcfs():
    args = parse_vcf_sources()
    annotate_source(
        source_path=args.source_path,
        output_name=args.output_name,
        cache_directory=args.cache_dir,
        assembly=args.assembly,
    )


if __name__ == "__main__":
    run_cli()
