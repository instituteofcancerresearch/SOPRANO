# Immunopeptidome builder

Users can build a custom immunopeptidome file based on a selection of HLA
alleles. This can be further refined by restricting or excluding Ensembl
transcript IDs.

```shell
soprano-hla2ip [-h] --alleles HLA_VALUES [HLA_VALUES ...] --output OUTPUT_ID [--cache CACHE_DIR] [--restrict [RESTRICTED_TRANSCRIPT_IDS ...]] [--excluded [EXCLUDED_TRANSCRIPT_IDS ...]]
```

`-h | --help`

:   Get help message

`-a | --alleles`

:   1 to 6 space separated HLA alleles, e.g., `HLA-A0201 HLA-A0101 HLA-C0701`

`-o | --output`

:   Name of the output immunopeptidome bed file. No extension required.

`-c | --cache`

:   Cache location of the output file. Defaults to
`app_sources/immunopeptidomes`

`-r | --restricted`

:   Space separated list of Ensembl transcript IDs to limit immunopeptidome to.

`-e | --excluded`

:   Space separated list of Ensembl transcript IDs to exclude from
immunopeptidome.
