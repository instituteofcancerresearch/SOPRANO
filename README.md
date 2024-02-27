[![SOPRANO (Dev) Tests](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/dev_tests.yml/badge.svg)](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/dev_tests.yml)

# SOPRANO: Selection On PRotein ANnotated regiOns
The SOPRANO method was developed to quantify selection in specific regions of 
the genome[^1]. 
It calculates ON- and OFF-target dN/dS 
using a set of annotated somatic point mutations and a transcript coordinates 
file. 

This repository is a Python reimplementation of the original method with two
main offerings:
1) An application interface to easily pre-process the relevant input files; and
2) A scalable command line utility to run the SOPRANO pipeline for dN/dS, 
serving local and high-performance computing environments, distributed over 
MPI.

[//]: # (_What inputs does SOPRANO require?_)

[//]: # ()
[//]: # (- A set of mutations &#40;missense/truncating&#41; and their respective functional )

[//]: # (annotation &#40;derived from ensembl VEP&#41;. Mutations can be filtered a priori by )

[//]: # (the user &#40;i.e. by predicted binding affinity or by expression status&#41;.)

[//]: # ()
[//]: # (- A set of ensembl transcript coordinates where selection will be estimated.)

[//]: # ()
[//]: # (_How does SOPRANO calculate dN/dS?_)

[//]: # ()
[//]: # (- "ON" dN/dS is the value calculated inside the coordinates provided using a )

[//]: # (192-trinucleotide correction signature obtained "on-the-fly" from the input )

[//]: # (mutation file. Alternatively, the user can provide a pre-calculated )

[//]: # (trinucleotide mutation frequency file. Importantly, ON dN/dS and OFF dN/dS )

[//]: # (&#40;the portion outside the coordinates provided&#41; will be calculated only using )

[//]: # (transcripts defined in this file. )

[//]: # ()
[//]: # (SOPRANO can be run from the command line, or via a streamlit app interface. )

[//]: # (There are additional tools built within the SOPRANO installation that enable )

[//]: # (users to download genomes, link existing VEP caches, and annotate VCF files.)

## Documentation

- [Installation](docs/INSTALL.md) - Building the environment and installing SOPRANO
- [CLI](docs/CLI.md) - Running SOPRANO from the command line.
- [GUI](docs/APP.md) - Running SOPRANO from your web browser.
- [Genomes](docs/GENOMES.md) - Downloading and preparing reference genome files.
- [Link VEP](docs/VEP.md) - Linking your existing VEP cache files.
- [Inputs 1](docs/INPUT.md) - Overview of the VEP annotated input.
- [Inputs 2](docs/BED.md) - Overview of the immunopeptidome input.
- [Output](docs/OUTPUT.md) - Overview of the summary statistics computed by SOPRANO.
- [Patient Specific](docs/PATIENTS.md) - How to obtain a patient specific immune dN/dS.
- [Notes](docs/NOTES.md) - Notes on and limitations of this software.

### Contact

Please raise issues and questions
[here](https://github.com/instituteofcancerresearch/SOPRANO/issues).

[^1]: See [Zapata et al](https://www.researchgate.net/publication/369116811_Immune_selection_determines_tumor_antigenicity_and_influences_response_to_checkpoint_inhibitors),
and the corresponding [original repository](https://github.com/luisgls/SOPRANO).