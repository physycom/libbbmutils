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
<a href="https://travis-ci.com/physycom/math_lib"> 
<div class="image">
<img src="https://travis-ci.com/physycom/math_lib.svg?token=ujjUseBa9hYbKckXBkxJ&branch=master" width="90" height="20" alt="Build Status"> 
</div>

## Purpose
This document presents a plain C library designed to handle vector and matrix algebra and 2x2 eigenvalue problems for symmetric matrices.
The library has been extensively used for our projects and is actively deployed into MetaSystem black boxes equipped with our algorithms to do some matrix calculations.

## Installation
Usually this library is deployed as a git submodule for other projects. It contains also some executables to do some self-tests of the routines.   
**make** and a **C99** compatible compiler are required for best results. Clone the repo and type ``make``, it should be enough in most cases!   
