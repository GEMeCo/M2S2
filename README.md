# M2S2

M2S2 is a simple C++ header-only library for linear algebra applied to mechanics of solids and structures. It is under development in Structural Engineering Department, at Sao Carlos School of Engineering, University of Sao Paulo.

<h1 align="center">
  <img alt="Banner" title="#Banner" height="150" src="./images/Icon_150.png" />
</h1>

## Copyright Information:
:ballot_box_with_check: This program is free software: you can redistribute it and/or modify it under the terms of the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

:ballot_box_with_check: Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer.

:ballot_box_with_check: Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer in the documentation and/or other materials provided with the distribution.

:ballot_box_with_check: This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
It is provided "AS IS". In no event shall the authors be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwire, arising from, out of or in connection with the software or the use of other dealing in the software.

:ballot_box_with_check: Neither the name of the copyright holder nor the names of any other contributors may be used to endorse or promote products derived from this software without specific prior written permission.

:triangular_ruler: None of the developers are **Software Engineers** of any sort. We are **Civil Engineers**, seeking to solve engineering problems. We know that this software has bugs and stuff to improve. We just don't have time to solve them all.

Documentation of this work, generated by doxygen, may also be used, distributed and modified, but under [Creative Commons Attribution-NonCommerical 4.0 International license](https://creativecommons.org/licenses/by-nc/4.0/).

<h1 align="center">
  <img alt="Banner" title="#Banner" height="100" src="./images/CC-BY-NC.jpg" />
</h1>

## :warning: Under development - Use at your own risk. Improper use will crash your application.

Table of Contents
=================
<!--ts-->
   * [About](#about)
   * [Citation](#citation)
   * [Features](#features)
   * [Resources](#resources)
   * [Building and Running](#how-to-run)
   * [Documentation](#documentation)
   * [Project Managers](#project)
   * [How to contribute](#how-to-contribute)
   * [Acknowledgement](#acknowledgement)
<!--te-->

About
-----
M2S2 is a simple C++ header only library for linear algebra applied to mechanics of solids and structures. It was developed primarly to help FEM coding.

Citation
--------
  Whether it was used in whole or parts, citation is a must!

 This library is in Zenodo - thus, we got a DOI: [doi.org/10.5281/zenodo.](https://doi.org/10.5281/zenodo.)
 
 (:construction: We'll get there, eventually :construction:).

Features
--------
- v0.1.0: Dyadics and Matrices
    - :boom: M2S2::Dyadic2S - Symmetric 2nd rank tensors of 2 or 3 dimensional vector space;
    - :boom: M2S2::Dyadic2N - Asymmetric 2nd rank tensors of 2 or 3 dimensional vector space;
    - :boom: M2S2::Dyadic4S - Symmetric 4th rank tensors of 2 or 3 dimensional vector space;
    - :boom: M2S2::Dyadic4C - 4th order tensor for orthotropic constitutive matrices, in either 2 or 3 dimensional vector space;
    - :boom: M2S2::MatrixS - Symmetric square matrix of any order.
    - :boom: M2S2::MatrixS - 2nd rank matrices of any size, even not square.

Resources
---------
M2S2 doesn't have any dependencies other than the C++ standard library.

Building and Running
--------------------
1. Clone the source code
`git clone --recursive https://github.com/GEMeCo/O2P2.git`

2. Just include the main header file (don't forget to setup including directories). That's all there is to it. Good luck :innocent:.
`#include "m2s2.h"`

Documentation
-------------
Software documentation is made directly from annotated sources by doxygen. Use doxywizard to create the documentation.

Project Managers
----------------
[Dorival Piedade Neto](http://lattes.cnpq.br/6930392733648456)

[Rodrigo Ribeiro Paccola](https://orcid.org/0000-0002-9228-4180)

[Rog�rio Carrazedo](https://orcid.org/0000-0003-2750-034X)


How to contribute
-----------------
```bash
Sorry, but we are not accepting external contributions, at least yet.
Nevertheless, you may use it as seen fit.
```

Acknowledgement
---------------
<h1 align="center">
  <img alt="Banner" title="#Banner" height="200" src="./images/logo_inst.png" />
</h1>

