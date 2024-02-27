# Installation

SOPRANO is compatible with Linux and MacOS. It is primarily written in Python,
but supplemented by some R, Perl, Bash, and bioinformatic utilities; such as
[bedtools](https://bedtools.readthedocs.io/en/latest/).

Due to the variety of languages, **SOPRANO is installed inside a conda
environment**. Note that this requires an existing Anaconda download, or
alternative variant, such as
[Mamba](https://mamba.readthedocs.io/en/latest/index.html); a highly efficient
re-implementation.

## Automatic installation

Users with **Anaconda**, or **Mamba**, available on their system can attempt an
automated installation from their command line:

```shell
. setup.sh <flavour>
```

The flavour of the installation is **optional**, and may be omitted. Currently
available options are:

- `mpi` - recommended for CLI users who user and HPC installations.
- `dev` - for developers; who are further recommended to `pre-commit install`
  hooks for standardized static analysis before CI.

Some users have reported issues with the automatic installation, in particular
on macOS systems. If you experience problems, please refer to the manual
installation. Pull requests in this area are highly welcome.

## Manual installation

This installation guide assumes that users are in the repository root.

1. Create and activate a host conda environment  
   ```shell
   conda create --name soprano
   conda activate soprano
   ```  

    1. Mac  
       ```shell
       conda config --env --set subdir osx-64
       conda env update ---file env.yml
       ```  

    2. Linux  
       ```shell
       conda env update --file env.yml
       ```  

2. Install R GitHub package dependencies  
   ```shell
   Rscript install_R_pkgs.R
   ```  

3. Install Python package and corresponding dependencies  
   ```shell
   pip install -e .[<flavour>]
   ```  

4. Decompress transcript files  
   ```shell
   gunzip -k data/aux_soprano/*.gz
   ```

## Validate your installation

After completing either installation procedure, users are encouraged to run
the installation tests, via

```shell
pytest -s tests/test_installation.py
```  

All tests should pass.
