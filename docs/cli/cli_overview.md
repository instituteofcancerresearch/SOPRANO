# CLI Overview

The SOPRANO CLI contains utilities for pipeline execution, in addition to
several utilities for preparing the necessary data inputs files.

## Core executables

- `soprano-run`  
Execute the SOPRANO pipeline to compute dN/dS, based on a choice of input
files and parameter choices.

- `soprano-app`  
Deploy the SOPRANO user interface to your local browser: run utilities and a
light-weight version of the pipeline.

## Additional utilities

- `soprano-link-vep`  
Link a VEP cache containing genome reference files to the SOPRANO data folders.

- `soprano-annotate`  
Annotate VCF files that can be handled by SOPRANO, describing somatic
mutations.

- `soprano-get-genome`  
Download an ensembl VEP genome reference file from the ensembl FTP server.

- `soprano-hla2ip`  
Prepare an immunopeptidome file based a discrete choice of HLA; and optional
transcript constraints.
