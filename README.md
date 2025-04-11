# README for the sdhdf\_tools repository

## About

This repository contains a collection of tools for working with 
Spectral-Domain Hierarchical Data Format (SDHDF) files.

**INSPECTA** is an 'INtegrated SDHDF Processing Engine in C for Telescope data Analysis'. 

Formerly known as 'sdhdfProc', INSPECTA is a software package to read, 
manipulate and process radio astronomy data in SDHDF format.

Please note: if you make use of this software, please reference Toomey et al. 2024:
https://doi.org/10.1016/j.ascom.2024.100804

**Author:**    George.Hobbs@csiro.au

**Copyright:** CSIRO 2020, 2021, 2022, 2023, 2024, 2025

## Contents

* containers
* docs
* inspecta 
* python (including pyINSPECTA, developed by Alec.Thompson@csiro.au)

## Building INSPECTA

To build INSPECTA, please follow the steps below:

1. Install prerequisites:

Erfa (https://github.com/liberfa/erfa) or the SOFA library (http://www.iausofa.org/)

HDF5 library (https://www.hdfgroup.org/downloads/hdf5) 

Pgplot library (http://www.astro.caltech.edu/~tjp/pgplot/) 

Calceph library (https://www.imcce.fr/recherche/equipes/asd/calceph/)

Cfitsio library

2. Navigate to the INSPECTA source code directory (inspecta/)

3. Run the bootstrap script:
```
./bootstrap
```

4. Run configure with e.g.:
```
./configure --prefix=/path/to/install CFLAGS=-I/path/to/include LDFLAGS=-L/path/to/lib
```
Where /path/to/install /path/to/include /path/to/lib are full paths to where you
wish to install it, full path to the prerequisite header files, 
and full path to the prerequisite libraries respectively.

5. Compile it with:
```
make
```

6. Install it with:
```
make install
```

7. Configure your SDHDF\_RUNTIME environmental variable to point to:

inspecta/runtime


## Branches

The repository has the following branches available:
* master
* sdhdfv3
Compatible with data written in SDHDF <=v3.0


## Test data

Please contact lawrence.toomey@csiro.au if you would like access to some tests data.


## More information

If you have any comments/queries, contact lawrence.toomey@csiro.au or george.hobbs@csiro.au

