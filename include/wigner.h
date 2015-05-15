/// \file wigner.h
/// \brief Header file for storing the Wigner state in Pauli representation
//
//	Copyright (c) 2014, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	References:
//	- Jianfeng Lu, Christian B. Mendl
//	  Numerical scheme for a spatially inhomogeneous matrix-valued quantum Boltzmann equation
//	  Journal of Computational Physics 291, 303-316 (2015)
//	  (arXiv:1408.1782)
//
//	- Martin L.R. F"urst, Christian B. Mendl, Herbert Spohn
//	  Matrix-valued Boltzmann equation for the Hubbard chain
//	  Physical Review E 86, 031122 (2012)
//	  (arXiv:1207.6926)
//_______________________________________________________________________________________________________________________
//

#ifndef WIGNER_H
#define WIGNER_H

#include "grid.h"


//_______________________________________________________________________________________________________________________
///
/// \brief Wigner state in physical velocity space discretized as [0:N/2-1,-N/2:-1]*(2L)/N in each direction
///
typedef struct
{
	grid_t comp[4];		//!< components in Pauli sigma representation
}
wignerV_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Wigner state in Fourier domain
///
typedef struct
{
	cgrid_t comp[4];	//!< components in Pauli sigma representation after Fourier transformation
}
wignerF_t;



#endif
