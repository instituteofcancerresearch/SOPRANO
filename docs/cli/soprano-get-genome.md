# Genome reference downloader

SOPRANO has a built-in genome reference downloader, that retrieves data
directly from the Ensembl ftp servers. This means that SOPRANO users do
not need to configure VEP, for example.

Users can execute the utility via

```shell
soprano-get-genome [-h] [--species SPECIES] [--assembly ASSEMBLY] [--release RELEASE] [--primary_assembly] [--download_only]
```

`-h | --help`

:   Get info on how to use this executable

`-s | --species`

:   Ensembl species Latin name. (default: homo_sapiens)

`-a | --assembly`

:   Ensembl genome assembly ID. (default: GRCh38)

`-r | --release`

:   Ensembl release number. (default: 110)

`-d | --download_only`

:   Don't decompress and compute fasta index file.

Downloads will be cached into the folder `ensembl_downloads/SPECIES/RELEASE`.
This mimics the VEP downloads structure. This is deliberate: users who
readily have genome references downloaded on their local laptop can leverage
the `soprano-link-vep` utility to crete a soft link to existing data in this
directory, which can then be used by SOPRANO.
