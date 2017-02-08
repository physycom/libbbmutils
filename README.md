---
documentclass: physycomen
title:  "math_lib"
author: "Fabbri, Sinigardi"
---

<a href="http://www.physycom.unibo.it"> 
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90" alt="© Physics of Complex Systems Laboratory - Physics and Astronomy Department - University of Bologna"> 
</div>
</a>

### Purpose
An **ANSI C** library to handle with vectors algebra, matrix algebra and eigenvalue problems. Its capabilities include:
- 2D, 3D, 6D vectors and their basic algebras, including rotations;
- 2D, 3D, 6D square matrices and their basic algebras, including rotations;
- 2D eigenvalue problems for symmetric matrices where exact diagonalization is possible.
The whole library is designed by making use of pointers in order to minimize the duplication of data with the aim of being employed in embedded systems. 
The library is also C++ compatible and comes equipped with a very basic test unit, try `make test`.

### Installation
**make** and a **C/C++** compatible compiler are required. Clone the repo and type ``make``, it should be enough in most cases to build the lib. To run tests type ``make test``.   
There's also a **VS2015** solution avalaible.   

### Documentation
The documentation can be generated by typing (``doxygen`` required) ``make doc``. To produce the pdf version you can uncomment the last two lines of the ``doc`` recipe in makefile.


