/// \file quadrature.h
/// \brief Header file for constructing the quadrature formulas required as input for the I1, I2, I3 and I4 integrals
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

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "grid.h"


//_______________________________________________________________________________________________________________________
///
/// \brief I1 quadrature formula representing Fourier transform of delta function or principal value
///
typedef struct
{
	unsigned int num;	//!< number of terms in sum
	grid2_t *psiR1;		//!< psiR1 weights
	grid_t  *psiR2;		//!< psiR2 weights
}
quadI1_t;


void QuadI1_Allocate(const unsigned int num, quadI1_t *quad);

void QuadI1_Delete(quadI1_t *quad);


//_______________________________________________________________________________________________________________________
///
/// \brief I2 quadrature formula representing Fourier transform of delta function
///
typedef struct
{
	unsigned int num;	//!< number of terms in sum
	grid_t *psiR1;		//!< psiR1 weights
	grid_t *psiR2;		//!< psiR2 weights
}
quadI2_t;


void QuadI2_Allocate(const unsigned int J, quadI2_t *quad);

void QuadI2_Delete(quadI2_t *quad);


//_______________________________________________________________________________________________________________________
///
/// \brief I3 quadrature formula representing Fourier transform of delta function
///
typedef struct
{
	unsigned int num;	//!< J, number of terms in sum
	cgrid_t *psiR1;		//!< psiR1 weights
	grid2_t *psiR2;		//!< psiR2 weights
}
quadI3_t;


void QuadI3_Allocate(const unsigned int num, quadI3_t *quad);

void QuadI3_Delete(quadI3_t *quad);


//_______________________________________________________________________________________________________________________
//


void FourierI1(const unsigned int J, const double L, const double R, quadI1_t *quad);

void FourierI2(const unsigned int J, const double L, const double R, quadI2_t *quad);

void FourierI3(const unsigned int J, const double L, const double R, const unsigned int M, quadI3_t *quad);

void FourierI4(const unsigned int J, const double L, const double R, quadI1_t *quad);



#endif
