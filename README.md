Simulation of a spatially inhomogeneous matrix-valued quantum Boltzmann equation using Fourier transformation
=============================================================================================================

C source code, test files and visualization routines for the simulations in "Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation" (see References).

Compiling the source code:
- Windows: a Visual Studio project file is provided in the *vcproj* folder
- Linux and similar: see the makefile in the *bin* folder; you may have to adapt paths according to your local installation

The *visualization* subfolder contains Matlab routines for visualizing the binary output files of the simulation code. Call `loadHomBinary()` and `loadInhomBinary('../bin/params_periodic_example.txt')` for visualizing the spatially homogeneous and inhomogeneous simulation examples, respectively.

License
-------
Copyright (c) 2014, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License
http://www.opensource.org/licenses/bsd-license.php


References
----------
1. Jianfeng Lu, Christian B. Mendl  
   Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation  
   Journal of Computational Physics 291, 303-316 (2015), [arXiv:1408.1782](http://arxiv.org/abs/1408.1782)
2. Martin L.R. Fürst, Christian B. Mendl, Herbert Spohn  
   Matrix-valued Boltzmann equation for the Hubbard chain  
   Physical Review E 86, 031122 (2012), [arXiv:1207.6926](http://arxiv.org/abs/1207.6926)
