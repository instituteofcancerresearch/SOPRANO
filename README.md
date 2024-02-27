[![SOPRANO (Dev) Tests](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/main_tests.yml/badge.svg)](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/main_tests.yml)
[![SOPRANO (Dev) Tests](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/dev_tests.yml/badge.svg)](https://github.com/instituteofcancerresearch/SOPRANO/actions/workflows/dev_tests.yml)

# SOPRANO: Selection On PRotein ANnotated regiOns
The SOPRANO method was developed to quantify selection in specific regions of 
the genome. 
It calculates ON- and OFF-target dN/dS 
using a set of annotated somatic point mutations and a transcript coordinates 
file. 

This repository is a Python reimplementation of the original method with two
main offerings:
1) An application interface to easily pre-process the relevant input files; and
2) A scalable command line utility to run the SOPRANO pipeline for dN/dS, 
serving local and high-performance computing environments, distributed over 
MPI.

### Documentation

Available on this [site](https://instituteofcancerresearch.github.io/SOPRANO/).

### Issues

Please raise issues and questions
[here](https://github.com/instituteofcancerresearch/SOPRANO/issues).

### References

See [Zapata et al](https://www.researchgate.net/publication/369116811_Immune_selection_determines_tumor_antigenicity_and_influences_response_to_checkpoint_inhibitors),
and the corresponding [original repository](https://github.com/luisgls/SOPRANO).