/// \file grid.h
/// \brief Header file for storing a square grid with real or complex entries in two dimensions
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

#ifndef GRID_H
#define GRID_H

#include <fftw3.h>

//_______________________________________________________________________________________________________________________
//


#define N_GRID		32				//!< number of grid points in each direction
#define N_MASK		(N_GRID - 1)	//!< bit mask for fast modulo N operation; N must be a power of 2


//_______________________________________________________________________________________________________________________
///
/// \brief Real grid defined as N x N matrix
///
typedef struct
{
	double data[N_GRID*N_GRID];				//!< data entries
}
grid_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Real grid defined as 2*N x 2*N matrix
///
typedef struct
{
	double data[4*N_GRID*N_GRID];			//!< data entries
}
grid2_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Complex grid defined as N x N matrix
///
typedef struct
{
	fftw_complex data[N_GRID*N_GRID];		//!< data entries
}
cgrid_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Complex grid defined as 2*N x 2*N matrix
///
typedef struct
{
	fftw_complex data[4*N_GRID*N_GRID];		//!< data entries
}
cgrid2_t;


//_______________________________________________________________________________________________________________________
//

void MultiplyPointwiseRR (const   grid_t *restrict_ f, const   grid_t *restrict_ g,   grid_t *restrict_ ret);
void MultiplyPointwiseRC (const   grid_t *restrict_ f, const  cgrid_t *restrict_ g,  cgrid_t *restrict_ ret);
void MultiplyPointwiseCC (const  cgrid_t *restrict_ f, const  cgrid_t *restrict_ g,  cgrid_t *restrict_ ret);
void MultiplyPointwiseRR2(const  grid2_t *restrict_ f, const  grid2_t *restrict_ g,  grid2_t *restrict_ ret);
void MultiplyPointwiseCC2(const cgrid2_t *restrict_ f, const cgrid2_t *restrict_ g, cgrid2_t *restrict_ ret);

void MultiplyPointwiseAssignRR (const   grid_t *restrict_ f,   grid_t *restrict_ ret);
void MultiplyPointwiseAssignRC (const   grid_t *restrict_ f,  cgrid_t *restrict_ ret);
void MultiplyPointwiseAssignRC2(const  grid2_t *restrict_ f, cgrid2_t *restrict_ ret);
void MultiplyPointwiseAssignCC2(const cgrid2_t *restrict_ f, cgrid2_t *restrict_ ret);


//_______________________________________________________________________________________________________________________
//

void GridZeroPad(const cgrid_t *restrict_ f, cgrid2_t *restrict_ ret);

void GridExtract(const cgrid2_t *restrict_ f, cgrid_t *restrict_ ret);



#endif
