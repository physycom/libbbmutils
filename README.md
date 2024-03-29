# libbbmutils

<a href="http://www.physycom.unibo.it">
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90" alt="© Physics of Complex Systems Laboratory - Physics and Astronomy Department - University of Bologna">
</div>
</a>

[![Physycom Continuous Integration](https://github.com/physycom/libbbmutils/actions/workflows/ccpp.yml/badge.svg)](https://github.com/physycom/libbbmutils/actions/workflows/ccpp.yml)

## Purpose

This document presents a C89 library designed to handle vector and matrix algebra and 2x2 eigenvalue problems for symmetric matrices.
The library has been extensively used for our projects and is actively deployed into MetaSystem black boxes equipped with our algorithms to do some matrix calculations.

## Details

A C89 library designed to handle problems with vectors algebra, matrix algebra and eigenvalues. Its capabilities include:

- 2D, 3D, 6D vectors and their basic algebras, including rotations;
- 2D, 3D, 6D square matrices and their basic algebras, including rotations;
- 2D eigenvalue problems for symmetric matrices where exact diagonalization is possible.

The whole library is designed by making use of pointers in order to minimize the duplication of data with the aim of being employed in embedded systems.
The library is also C++ compatible and comes equipped with a very basic test unit. After building executables, run them with

```pwsh
./ci/verify-test-log.ps1 -TestBinFolder ./bin/ -TestLogFolder ./test/
```

## Installation

**CMake**, a **C89** and a **C++11** compatible compiler are required. To build the executable, clone the repo and then type  

```pwsh
./ci/build.ps1 -UseVCPKG
```

If you want to build also documentation, please run

```pwsh
./ci/build.ps1 -UseVCPKG -BuildDocumentation
```
