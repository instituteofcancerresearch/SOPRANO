# Mutation annotations

SOPRANO has a specific format which is required to describe somatic mutations.

The format is most easily derived from a VCF file, or collection of VCF files.
Users should use

```shell
soprano-annotate [-h] [--source SOURCE_PATH] [--output OUTPUT_NAME] [--cache CACHE_DIR] [--assembly ASSEMBLY]
```

`-h | --help`

:   Get help message.

`-s | --source`

:   Provide the path to a single VCF file, or a directory containing VCF files.

`-o | --output`

:   Provide a name for the output file. If not provided, when applied to as
single input file, will default to the input file name with a `.anno`
extension. If multiple files are selected (via prescribing a directory), this
will field is required.

`-d | --cache`

:   Provide a path to a directory to cache the output inside of. Defaults to
`/app_sources/annotated_inputs`.

`-a | --assembly`

:   Provide a genome assembly that is compatible with the input VCF source
files. By default, this assumed to be GRCh38.
