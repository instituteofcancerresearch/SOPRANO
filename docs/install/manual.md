# Manual installation

This installation guide assumes that users are in the repository root.

1. **Create and activate a host conda environment**

        conda create --name soprano
        conda activate soprano

    1. **macOS**

            conda config --env --set subdir osx-64
            conda env update ---file env.yml

    2. **Linux**

            conda env update --file env.yml

2. **Install R GitHub package dependencies**

        Rscript install_R_pkgs.R

3. **Install Python package and corresponding dependencies**

        pip install -e .[<flavour>]

   See the previous page for flavor options; omit braces `[...]` if not
   required.

4. **Decompress transcript files**

        ginzip -k data/aux_soprano/*.gz

5. **Validate your installation**

        pytest -s tests/test_installation.py

All tests should pass at the end of this installation.

[//]: # (### Conda environment)

[//]: # ()

[//]: # (The most straightforward way to manage these dependencies is to create a)

[//]: # (`conda` environment. From the repository root:)

[//]: # ()

[//]: # (```shell)

[//]: # (conda env create -f env.yml)

[//]: # (```)

[//]: # ()

[//]: # (After creation, the environment can be activated with)

[//]: # ()

[//]: # (```shell)

[//]: # (conda activate soprano)

[//]: # (```)

[//]: # ()

[//]: # ([//]: # &#40;#### MacOS note&#41;)

[//]: # ()

[//]: # (Mac users with the latest system architecture may instead be required to:)

[//]: # ()

[//]: # (```shell)

[//]: # (conda create --name soprano)

[//]: # (conda activate soprano)

[//]: # (conda config --env --set subdir osx-64)

[//]: # (conda env update --file $DEPS_FILE)

[//]: # (```)

[//]: # ()

[//]: # (```shell)

[//]: # (conda config --env --set subdir osx-64)

[//]: # (```)

[//]: # ()

[//]: # (### Non-conda users)

[//]: # ()

[//]: # (Users who do not wish to use `conda` should ensure that they have appropriate)

[//]: # (versions of `Python`, `R`, `Perl` and `Bedtools` installed on their system.)

[//]: # ()

[//]: # (To install SOPRANO as a Python package, it is recommended to create a Python)

[//]: # (virtual environment, e.g., from the repository root)

[//]: # ()

[//]: # (```shell)

[//]: # (python3 -m venv .venv)

[//]: # (```)

[//]: # ()

[//]: # (After creation, the environment can be activated with)

[//]: # ()

[//]: # (```shell)

[//]: # (source .venv/bin/activate)

[//]: # (```)

[//]: # ()

[//]: # (#### MacOS note)

[//]: # ()

[//]: # (Some GNU command line utilities are not shipped with MacOS natively, but can)

[//]: # (can be installed with e.g.,)

[//]: # ()

[//]: # (```shell)

[//]: # (brew install coreutils)

[//]: # (```)

[//]: # ()

[//]: # (### Installing R package dependencies)

[//]: # ()

[//]: # (To install the dependencies for mutation annotations, from the repository)

[//]: # (root)

[//]: # ()

[//]: # (```shell)

[//]: # (Rscript install_R_pkgs.R)

[//]: # (```)

[//]: # ()

[//]: # (Conda users should readily have most of these packages installed from the)

[//]: # (previous environment build, though `dndscv` will still be installed directly)

[//]: # (from GitHub.)

[//]: # ()

[//]: # (### Installing the SOPRANO Python package)

[//]: # ()

[//]: # (Ensure an appropriate environment is activated, then)

[//]: # ()

[//]: # (```shell)

[//]: # (pip install -e .)

[//]: # (```)

[//]: # ()

[//]: # (This will install a standard &#40;though editable version of SOPRANO&#41;.)

[//]: # ()

[//]: # (### Decompressing transcripts)

[//]: # ()

[//]: # (Ensure that transcript data is decompressed via)

[//]: # ()

[//]: # (```shell)

[//]: # (gzip -k "src/SOPRANO/data/ensemble_transcriptID.fasta.gz")

[//]: # (```)

[//]: # ()

[//]: # (### Validating the installation)

[//]: # ()

[//]: # (To check that the validity of your installation, users can run)

[//]: # ()

[//]: # (```shell)

[//]: # (pytest -s tests/test_installation.py)

[//]: # (```)

[//]: # ()

[//]: # (### Developer notes)

[//]: # ()

[//]: # (Developers should install the SOPRANO package with additional dependencies)

[//]: # ()

[//]: # (```shell)

[//]: # (pip install -e .[dev])

[//]: # (```)

[//]: # ()

[//]: # (This will install additional packages &#40;such as black, ruff, isort and pytest&#41;.)

[//]: # (Consistent styling can be enforced with pre-commit hooks: `pre-commit install`.)
