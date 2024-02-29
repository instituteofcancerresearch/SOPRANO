# SOPRANO App

SOPRANO has a browser-based user interface powered by `streamlit`. This can
be deployed from the command line via `soprano-app`. This will present the
user with a visual way of interacting with the CLI utilities. It also
offers a light-weight version of the pipeline; though this is restricted to
single-sample pipeline runs, with no KDE as described in the CLI.

## Data sources

By default, the application will serve options for the example annotated
mutation and immunopeptidome files shipped with SOPRANO. Additional files
can be detected by the application by their placement in the `./app_sources`
folders:

`./app_sources/annotated_inputs`

:   VEP annotated mutation files placed in this directory will be detected,
so long as they have the extension pattern `*anno*`. E.g., `mutations.anno`
or `mutations.annotated`.

`./app_sources/immunopeptidomes`

:   User defined immunopeptidomes BED files placed in this directory will be
detected, so long as they have the extension `.bed`. E.g., `immuno.bed`.

`./app_sources/coordinate_files`

:   User defined BED files that can be used for randomization will be detected,
so long as they have the extension `.bed`. E.g., `randoms.bed`.

**Note:** Quantities generated on-the-fly during an app session will be cached
into these folders.

***

## Serving the application via Docker

The application can also be constructed as a docker image, and served more
generally this way. This can be manually constructed from the repository root
via

```shell
docker build -t soprano .
```

To interact with data sources, users should bind mount whilst running the
application, for example, with

```shell
docker run -d -p 8501:8501 --name SOPRANO_APP -v "$(pwd)"/ensembl_downloads:/app/ensembl_downloads soprano
```

***

### Public content

Alternatively, there is a
[publicly available image](https://hub.docker.com/r/icrsc/soprano) that can
be pulled down directly, and executed in the same way. A version of this
application is running at the ICR that you can view
[here](https://software.icr.ac.uk/app/soprano).

### Notes on the containerized app

Data in general will not persist between sessions. For example, an annotated
mutation file will not be saved directly onto your file system. Users are
therefore encouraged to download any results that they want to use via the
appropriate buttons.
