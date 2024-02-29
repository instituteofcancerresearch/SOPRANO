# VEP cache

Genome reference files can consume large amounts of storage. Many users of
SOPRANO may have existing VEP caches with genome references therein, that could
be used directly by SOPRANO.

To this end, users can link their existing reference files to the SOPRANO
Ensembl downloads folder `ensembl_downloads/`.

Users should run:

```shell
soprano-link-vep [-h] [--cache SRC_CACHE]
```

`-h | --help`

:   Get help output

`-c | --cache`

:   Provide the path to the ensembl vep cache. (default: /home/USER/.vep)

The cache argument is useful if you have cached VEP data into a non-standard
location.
