import json
import logging
import pathlib
import warnings
from argparse import Namespace
from dataclasses import dataclass
from typing import Set

import pandas as pd

from SOPRANO.utils.mpi_utils import (
    COMM,
    RANK,
    as_single_process,
    item_selection,
)
from SOPRANO.utils.path_utils import Directories, genome_pars_to_paths
from SOPRANO.utils.url_utils import (
    build_ensembl_urls,
    check_ensembl_file_url,
    compute_chrom_sizes,
    compute_fasta_index,
    decompress,
    download_from_url,
    filename_from_url,
    find_earliest_release,
    find_latest_release,
)


@dataclass(frozen=True)
class TranscriptPaths:
    transcript_length: pathlib.Path
    protein_transcript_length: pathlib.Path
    transcript_fasta: pathlib.Path

    @classmethod
    def defaults(cls):
        return cls(
            transcript_length=Directories.soprano_aux_files(
                "ensemble_transcript.length"
            ),
            protein_transcript_length=Directories.soprano_aux_files(
                "ensemble_transcript_protein.length"
            ),
            transcript_fasta=Directories.soprano_aux_files(
                "ensemble_transcriptID.fasta"
            ),
        )


@dataclass(frozen=True)
class GenomePaths:
    sizes: pathlib.Path
    fasta: pathlib.Path

    @classmethod
    def GRCh37(cls, release=110):
        fasta, sizes = genome_pars_to_paths("GRCh37", release)
        return cls(sizes=sizes, fasta=fasta)

    @classmethod
    def GRCh38(cls, release=110):
        fasta, sizes = genome_pars_to_paths("GRCh38", release)
        return cls(sizes=sizes, fasta=fasta)


@dataclass(frozen=True)
class AuxiliaryPaths:
    genes_to_exclude: pathlib.Path
    intron_length: pathlib.Path

    @classmethod
    def defaults(cls):
        return cls(
            genes_to_exclude=Directories.soprano_aux_files(
                "genes2exclude.txt"
            ),
            intron_length=Directories.soprano_aux_files(
                "transcript_intron_length.bed"
            ),
        )


class AnalysisPaths:
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        random_regions: pathlib.Path | None = None,
    ):
        self.analysis_name = analysis_name
        self.input_path = input_path
        self.bed_path = bed_path
        self.random_regions_path = random_regions
        self.cache_dir = cache_dir

        # Transcripts
        self.filtered_protein_transcript = self._cached_path(
            "protein_length_filt", "txt"
        )
        self.filtered_transcript = self._cached_path(
            "transcript_length_filt", "txt"
        )

        # Exclusion regions for randomization
        self.exclusions = self._cached_path("exclusion", "ori")
        self.exclusions_sorted = self._cached_path("exclusion", "bed")
        self.exclusions_shuffled = self._cached_path("epitopes", "ori2")

        # Epitope files
        self.epitopes = self._cached_path("epitopes", "bed")
        self.epitopes_cds = self._cached_path("epitopes_cds", "bed")
        self.epitopes_cds_fasta = self._cached_path("epitopes_cds", "fasta")
        self.epitopes_trans_regs = self._cached_path(
            "epitopes_cds", "transcript_regions"
        )
        self.epitopes_trans_regs_sum = self._cached_path(
            "epitopes_cds", "transcript_regions", "sum"
        )

        # Complement (intra) epitope files
        self.intra_epitopes = self._cached_path("intra_epitopes", "bed")
        self.intra_epitopes_tmp = self._cached_path(
            "intra_epitopes", "bed", "tmp"
        )
        self.intra_epitopes_prot = self._cached_path(
            "intra_epitopes_prot", "bed"
        )
        self.intra_epitopes_prot_tmp = self._cached_path(
            "intra_epitopes_prot", "bed", "tmp"
        )
        self.intra_epitopes_cds = self._cached_path("intra_epitopes_cds")
        self.intra_epitopes_cds_fasta = self._cached_path(
            "intra_epitopes_cds", "fasta"
        )
        self.intra_epitopes_trans_regs = self._cached_path(
            "intra_epitopes_cds", "transcript_regions"
        )
        self.intra_epitopes_trans_regs_sum = self._cached_path(
            "intra_epitopes_cds", "transcript_regions", "sum"
        )

        # Analysis files
        self.sim_fixed = self._cached_path("sim_fixed")
        self.col_corrected = self._cached_path("col_corrected")
        self.contextualised = self._cached_path("contextualised")
        self.flagged = self._cached_path("flagged")
        self.triplet_counts = self._cached_path("triplets", "counts")
        self.final_epitope_corrections = self._cached_path(
            "corrected_matrix", "epitopes"
        )
        self.final_intra_epitope_corrections = self._cached_path(
            "corrected_matrix", "intra_epitopes"
        )
        self.epitope_nans = self._cached_path("epitopes", "nans")
        self.intra_epitope_nans = self._cached_path("intra_epitopes", "nans")

        # variant counts
        self.variants_silent = self._cached_path("variants", "silent", "bed")
        self.variants_nonsilent = self._cached_path(
            "variants", "nonsilent", "bed"
        )
        self.variants_missense = self._cached_path(
            "variants", "missense", "bed"
        )
        self.variants_intronic = self._cached_path(
            "variants", "intronic", "bed"
        )
        self.raw_silent_count = self._cached_path("raw", "silent", "count")
        self.raw_nonsilent_count = self._cached_path(
            "raw", "nonsilent", "count"
        )
        self.raw_missense_count = self._cached_path("raw", "missense", "count")
        self.in_silent_count = self._cached_path("in", "silent", "count")
        self.in_nonsilent_count = self._cached_path("in", "nonsilent", "count")
        self.in_missense_count = self._cached_path("in", "missense", "count")
        self.out_silent_count = self._cached_path("out", "silent", "count")
        self.out_nonsilent_count = self._cached_path(
            "out", "nonsilent", "count"
        )
        self.out_missense_count = self._cached_path("out", "missense", "count")

        self.data_epitopes = self._cached_path("data", "epitopes")

        self.intron_rate = self._cached_path("intron", "rate")

        self.results_path = self._cached_path("results", "tsv")

        self.log_path = self._cached_path("log")

    def _cached_path(self, *extensions):
        file_name = f"{self.analysis_name}"

        if len(extensions) > 0:
            file_name += "." + ".".join([*extensions])

        return self.cache_dir.joinpath(file_name)

    def is_complete(self):
        return self.results_path.exists()


_NAMESPACE_KEYS = (
    "analysis_name",
    "input_path",
    "bed_path",
    "cache_dir",
    "random_regions",
    "transcript",
    "protein_transcript",
    "transcript_ids",
    "use_ssb192",
    "keep_drivers",
    "seed",
    "species",
    "assembly",
    "release",
    "n_samples",
)


def check_cache_path(cache_dir: pathlib.Path, name: str) -> pathlib.Path:
    if cache_dir.exists() and cache_dir.is_dir():
        job_cache = cache_dir.joinpath(name)
        job_cache.mkdir(exist_ok=True)
        return job_cache
    else:
        raise NotADirectoryError(cache_dir)


def init_logger(name: str, log_path: pathlib.Path):
    log_level = logging.INFO
    proc_rank = "%04d" % RANK
    log_format = f"%(asctime)s | proc {proc_rank} | %(message)s"
    formatter = logging.Formatter(log_format)

    logger = logging.getLogger(name)
    logger.setLevel(log_level)

    log_file_handler = logging.FileHandler(log_path.as_posix())
    log_file_handler.setLevel(log_level)
    log_file_handler.setFormatter(formatter)

    logger.addHandler(log_file_handler)

    return logger


class Parameters(AnalysisPaths):
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        random_regions: pathlib.Path | None,
        use_ssb192: bool,
        use_random: bool,
        exclude_drivers: bool,
        seed: int | None,
        transcripts: TranscriptPaths,
        genomes: GenomePaths,
    ):
        super().__init__(
            analysis_name, input_path, bed_path, cache_dir, random_regions
        )

        self.transcripts = transcripts
        self.genomes = genomes
        self.use_ssb192 = use_ssb192
        self.use_random_regions = random_regions is not None
        self.use_random = use_random
        self.exclude_drivers = exclude_drivers
        self.seed = seed
        self.logger = init_logger(self.analysis_name, self.log_path)

        self.log("parameters initialized")

    def log(self, msg: str) -> None:
        self.logger.info(msg)


class GlobalParameters:
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        job_cache: pathlib.Path,
        random_regions: pathlib.Path | None,
        use_ssb192: bool,
        exclude_drivers: bool,
        seed: int | None,
        transcripts: TranscriptPaths,
        genomes: GenomePaths,
        n_samples: int,
    ):
        # Sanitized
        self.job_cache = check_cache_path(job_cache, analysis_name)
        self.seed = GlobalParameters.check_seed(seed)

        # Assumed OK
        self.input_path = input_path
        self.bed_path = bed_path
        self.random_regions = random_regions
        self.use_ssb192 = use_ssb192
        self.exclude_drivers = exclude_drivers
        self.transcripts = transcripts
        self.genomes = genomes
        self.n_samples = n_samples

        self.get_all_samples(_init=True)
        self.cache_ordered_params()

        # Post processing data files
        self.samples_path = self.job_cache.joinpath("samples_df.csv")
        self.samples_meta_path = self.job_cache.joinpath("samples_df.meta")

    def get_ordered_params(self):
        kwargs = self.__dict__.copy()

        def expand_kwargs(d: dict):
            _rerun = False

            for key, value in d.items():
                if hasattr(value, "__dict__"):
                    d2 = value.__dict__
                    del d[key]
                    d.update(d2)
                    _rerun = True
                    break

            if _rerun:
                return expand_kwargs(d)

            return d

        del kwargs["n_samples"]

        expanded = expand_kwargs(kwargs)

        return {k: str(kwargs[k]) for k in sorted(expanded)}

    def get_params_path(self):
        return self.job_cache / "pipeline.params"

    @as_single_process()
    def cache_ordered_params(self):
        params_path = self.get_params_path()

        if not params_path.exists():
            ordered_params = self.get_ordered_params()

            with open(params_path, "w") as f:
                json.dump(ordered_params, f, indent=4)
        else:
            self.check_params()

    def check_params(self):
        current_params = self.get_ordered_params()
        cached_params_path = self.get_params_path()

        with open(cached_params_path, "r") as f:
            cached_params = json.load(f)

        current_param_keys = tuple(current_params.keys())
        cached_param_keys = tuple(cached_params.keys())

        if current_param_keys != cached_param_keys:
            raise KeyError(f"{current_param_keys} != {cached_param_keys}")

        for _key in current_param_keys:
            current_value = current_params[_key]
            cached_value = cached_params[_key]

            if current_value != cached_value:
                raise ValueError(
                    f"Value mismatch @ {_key}: "
                    f"{current_value} != {cached_value}.\n"
                    f"Different pipeline runs must have distinct directories!"
                )

    @as_single_process()
    def gather(self):
        sample_results_paths = [
            self.get_sample(idx).results_path for idx in range(self.n_samples)
        ]

        for expected_results_path in sample_results_paths:
            if not expected_results_path.exists():
                warnings.warn(
                    f"Expected results file not found: "
                    f"{expected_results_path} ... \n"
                    f"omitting from post-processing."
                )

                sample_results_paths.remove(expected_results_path)

        if len(sample_results_paths) == 0:
            raise ValueError(f"No sample results found for {self.job_cache}.")

        joined_df: pd.DataFrame | None = None

        with open(self.samples_meta_path, "w") as f:
            for path in sample_results_paths:
                if joined_df is None:
                    joined_df = pd.read_csv(path, sep="\t")
                else:
                    joined_df = pd.concat(
                        [joined_df, pd.read_csv(path, sep="\t")],
                        ignore_index=True,
                    )

                f.write(f"{path.as_posix()}\n")

        # Dropped estimateed statistics... don't mean much in this context
        joined_df.drop(
            columns=[
                "ON_Low_CI",
                "ON_High_CI",
                "OFF_Low_CI",
                "OFF_High_CI",
                "Pvalue",
            ]
        )

        joined_df.to_csv(self.samples_path)

    @staticmethod
    def split_joined_df(
        joined_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        exonic_only = joined_df[["Coverage"] == "Exonic_Only"]
        exonic_intronic = joined_df[["Coverage"] == "Exonic_Intronic"]

        return exonic_only, exonic_intronic

    @staticmethod
    def check_seed(seed: int | None) -> int:
        if seed is None or seed < 0:
            print(f"WARNING: Random seed={seed}<0, assigning default: 1234.")
            return 1234

        return seed

    @classmethod
    def from_namespace(cls, namespace: Namespace):
        input_namespace_keys = namespace.__dict__.keys()
        if set(_NAMESPACE_KEYS) != set(input_namespace_keys):
            raise KeyError(
                f"{sorted(_NAMESPACE_KEYS)} != {sorted(input_namespace_keys)}"
            )

        transcripts = TranscriptPaths(
            namespace.transcript,
            namespace.protein_transcript,
            namespace.transcript_ids,
        )

        species = namespace.species
        assembly = namespace.assembly
        release = namespace.release

        if assembly == "GRCh37":
            assert species == "homo_sapiens", species
            genomes = (
                EnsemblData.homo_sapiens_GRCh37().get_genome_reference_paths(
                    release
                )
            )
        else:
            genomes = EnsemblData(
                species=species, assembly=assembly
            ).get_genome_reference_paths(release)

        n_samples = namespace.n_samples

        return cls(
            analysis_name=namespace.analysis_name,
            input_path=namespace.input_path,
            bed_path=namespace.bed_path,
            job_cache=namespace.cache_dir,
            random_regions=namespace.random_regions,
            use_ssb192=namespace.use_ssb192,
            exclude_drivers=not namespace.keep_drivers,
            seed=namespace.seed,
            transcripts=transcripts,
            genomes=genomes,
            n_samples=n_samples,
        )

    def get_data(self):
        return self.get_sample(-1)

    def get_sample(self, idx: int, _init=False):
        if not (-1 <= idx < self.n_samples):
            raise ValueError(
                f"Index {idx} is out of range for number of samples: "
                f"{self.n_samples}"
            )

        sample_seed = self.seed + idx if idx > -1 else None
        subdir_name = "data" if idx == -1 else "sample_%04d" % idx
        sample_cache = self.job_cache / subdir_name
        use_random = idx > -1

        if _init:
            if COMM.Get_rank() == 0:
                sample_cache.mkdir(parents=True, exist_ok=True)
            COMM.Barrier()
            return

        sample_kwargs = self.__dict__.copy()
        del sample_kwargs["n_samples"]
        del sample_kwargs["job_cache"]
        sample_kwargs["seed"] = sample_seed
        sample_kwargs["cache_dir"] = sample_cache
        sample_kwargs["analysis_name"] = subdir_name
        sample_kwargs["use_random"] = use_random

        return Parameters(**sample_kwargs)

    def get_all_samples(self, _init=False):
        samples = [
            self.get_sample(idx, _init=_init)
            for idx in range(-1, self.n_samples)
        ]

        if not _init:
            return [s for s in samples if not s.is_complete()]

    def get_worker_samples(self):
        worker_samples = item_selection(self.get_all_samples())
        return worker_samples


class SOPRANOError(Exception):
    pass


class _GatherReferences:
    # Urls
    toplevel_url: str
    primary_assembly_url: str

    # Params
    species: str
    assembly: str

    # Status
    toplevel_gz_done: Set[int] = set()
    toplevel_fa_done: Set[int] = set()
    toplevel_fai_done: Set[int] = set()
    primary_assembly_gz_done: Set[int] = set()
    primary_assembly_fa_done: Set[int] = set()
    primary_assembly_fai_done: Set[int] = set()
    sizes_done: Set[int] = set()

    def _dest_directory(self, release: int):
        return (
            Directories.ensembl_downloads(self.species)
            / f"{release}_{self.assembly}"
        )

    def _dest_fa_gz(self, release: int, _toplevel: bool):
        return self._dest_directory(release) / filename_from_url(
            self.toplevel_url if _toplevel else self.primary_assembly_url
        ).format(RELEASE=release)

    def _dest_fa(self, release: int, _toplevel: bool):
        return self._dest_fa_gz(release, _toplevel).with_suffix("")

    def _dest_fai(self, release: int, _toplevel: bool):
        return self._dest_fa_gz(release, _toplevel).with_suffix(".fai")

    def dest_chrom(self, release: int, _toplevel: bool):
        return self._dest_fa(release, _toplevel).with_suffix(".chrom")

    def toplevel_fa_gz_path(self, release: int):
        return self._dest_fa_gz(release, _toplevel=True)

    def toplevel_fa_path(self, release: int):
        return self._dest_fa(release, _toplevel=True)

    def toplevel_fai_path(self, release: int):
        return self._dest_fai(release, _toplevel=True)

    def toplevel_chrom_path(self, release: int):
        return self.toplevel_fa_path(release).with_suffix(".chrom")

    def primary_assembly_fa_gz_path(self, release: int):
        return self._dest_fa_gz(release, _toplevel=False)

    def primary_assembly_fa_path(self, release: int):
        return self._dest_fa(release, _toplevel=False)

    def primary_assembly_fai_path(self, release: int):
        return self._dest_fai(release, _toplevel=False)

    def _download(self, release: int, _toplevel):
        if _toplevel:
            source_url = self.toplevel_url
            dest_path = self.toplevel_fa_gz_path(release)
            decompressed_path = self.toplevel_fa_path(release)
        else:
            source_url = self.primary_assembly_url
            dest_path = self.primary_assembly_fa_gz_path(release)
            decompressed_path = self.primary_assembly_fa_path(release)

        if not (decompressed_path.exists() or dest_path.exists()):
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            check_ensembl_file_url(source_url, release)
            download_from_url(
                source_url.format(RELEASE=release),
                target_path=dest_path,
            )

    def _check_release_ok(self, release):
        min_release = find_earliest_release(self.toplevel_url)
        max_release = find_latest_release(self.toplevel_url)

        if not (min_release <= release <= max_release):
            raise ValueError(release)

    def download_toplevel(self, release):
        if release not in self.toplevel_gz_done:
            self._check_release_ok(release)

            if not self.toplevel_fa_gz_path(release).exists():
                self._download(release, _toplevel=True)

            self.toplevel_gz_done.add(release)

    def download_primary_assembly(self, release):
        if release not in self.primary_assembly_gz_done:
            self._check_release_ok(release)

            if not self.primary_assembly_fa_gz_path(release).exists():
                self._download(release, _toplevel=False)

            self.primary_assembly_gz_done.add(release)

    def decompress_toplevel(self, release):
        if release not in self.toplevel_fa_done:
            if not self.toplevel_fa_path(release).exists():
                decompress(self.toplevel_fa_gz_path(release))

            self.toplevel_fa_done.add(release)

    def decompress_primary_assembly(self, release):
        if release not in self.primary_assembly_fa_done:
            if not self.primary_assembly_fa_path(release).exists():
                decompress(self.primary_assembly_fa_gz_path(release))

            self.primary_assembly_fa_done.add(release)

    def compute_chrom_sizes(self, release):
        if release not in self.sizes_done:
            if not self.toplevel_chrom_path(release).exists():
                compute_chrom_sizes(self.toplevel_fai_path(release))

            self.sizes_done.add(release)

    def compute_fasta_index_toplevel(self, release):
        if release not in self.toplevel_fai_done:
            if not self.toplevel_fai_path(release).exists():
                compute_fasta_index(self.toplevel_fa_path(release))

            self.toplevel_fai_done.add(release)

    def compute_fasta_index_primary_assembly(self, release):
        if release not in self.primary_assembly_fai_done:
            if not self.primary_assembly_fai_path(release).exists():
                compute_fasta_index(self.primary_assembly_fa_path(release))

            self.primary_assembly_fai_done.add(release)

    def compute_all_toplevel(self, release):
        self.download_toplevel(release)
        self.decompress_toplevel(release)
        self.compute_fasta_index_toplevel(release)
        self.compute_chrom_sizes(release)

    def compute_all_primary_assembly(self, release):
        self.download_primary_assembly(release)
        self.decompress_primary_assembly(release)
        self.compute_fasta_index_primary_assembly(release)

    def get_genome_reference_paths(self, release):
        return GenomePaths(
            sizes=self.toplevel_chrom_path(release),
            fasta=self.toplevel_fa_path(release),
        )


class EnsemblData(_GatherReferences):
    def __init__(self, species: str, assembly: str, _init_urls=True):
        self.species = species
        self.assembly = assembly

        if _init_urls:
            url_dict = build_ensembl_urls(species, assembly)
            self.toplevel_url = url_dict["toplevel"]
            self.primary_assembly_url = url_dict["primary_assembly"]

    @classmethod
    def homo_sapiens_GRCh38(cls):
        return cls("homo_sapiens", "GRCh38")

    @classmethod
    def homo_sapiens_GRCh37(cls):
        # GRCh37 is actually has a deviant url structure, so manually set here
        toplevel_url = (
            "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
        )

        primary_assembly_url = (
            "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
        )

        species = "homo_sapiens"
        assembly = "GRCh37"

        obj = cls(species, assembly, _init_urls=False)
        obj.toplevel_url = toplevel_url
        obj.primary_assembly_url = primary_assembly_url
        return obj
