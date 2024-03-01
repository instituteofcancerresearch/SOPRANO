# SOPRANO

**_Selection On PRotein ANnotated regiOns_**

***

[![Main Tests](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/main_tests.yml/badge.svg?branch=python)](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/main_tests.yml)
[![Dev Tests](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/dev_tests.yml/badge.svg?branch=python-dev)](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/dev_tests.yml)
[![Python 3.10|11](https://img.shields.io/badge/Python-3.10_|_3.11-blue.svg)](https://img.shields.io/badge/Python-3.10_|_3.11-blue.svg)

The SOPRANO method was developed to quantify selection in specific regions of
the genome.
It calculates ON- and OFF-target dN/dS
using a set of annotated somatic point mutations and a transcript coordinates
file.

This repository is a Python reimplementation of the original method with two
main offerings:

1. An application interface to easily pre-process the relevant input files.  
2. A scalable command line utility to run the SOPRANO pipeline for dN/dS,
serving local and high-performance computing environments, distributed over
MPI.

## Quick links

- **Documentation** is available on this [site](https://instituteofcancerresearch.github.io/SOPRANO/).

- **Issues** should be raised
[here](https://github.com/instituteofcancerresearch/SOPRANO/issues).

- **Reference** the SOPRANO method with [Zapata et al](https://www.researchgate.net/publication/369116811_Immune_selection_determines_tumor_antigenicity_and_influences_response_to_checkpoint_inhibitors).

- **View** the original repository [here](https://github.com/luisgls/SOPRANO).
