# Input mutation files

This data corresponds to the SOPRANO CLI argument `-i | --input`. Via the
app, this appears under "Select a VEP annotated file".

## VEP annotation

Mutations should be annotated with Variant Effect Predictor (VEP) using the
[default output format](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html).

The input mutation file is the output from VEP, obtained via the
command

**_This command was deprecated_**

assuming the input is in the default Ensembl input format.

## Putative germline variants

If you want to filter putative germline variants, use the option
`--plugin ExAC`
when running VEP. It is important that you restrict your analysis to the list
of [Ensembl transcripts](../data/aux_soprano/ensemble_transcriptID.fasta)
observed in SOPRANO `src/SOPRANO/data` folder. Be aware of updates on the
EnsemblID from previous versions.

## Translating VCF annotated files

To convert a VCF annotated file to the input format for SOPRANO, the user can
run [vep-annotation-reporter](https://vatools.readthedocs.io/en/latest/vep_annotation_reporter.html)
using the command

**_This command was deprecated_**

## Examples

Example input files are bundled with SOPRANO:

- [TCGA-05-4396-01A-21D-1855-08.annotated](../data/example_annotations/TCGA-05-4396-01A-21D-1855-08.annotated)
- [rand100K.annotated](../data/example_annotations/rand100K.annotated)

More examples can be found on synapse: ID syn11681983.
