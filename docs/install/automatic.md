# Installation

SOPRANO is compatible with Linux and MacOS. It is primarily written in Python,
but supplemented by some R, Perl, Bash, and bioinformatic utilities; such as
[bedtools](https://bedtools.readthedocs.io/en/latest/).

Due to the variety of languages, **installation is recommended via a
conda environment**. Note that this requires an existing Anaconda download, or
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
