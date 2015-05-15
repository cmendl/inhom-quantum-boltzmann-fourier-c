/// \file grid.c
/// \brief Operations on square grids with real or complex entries in two dimensions
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

#include "grid.h"
#include <memory.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'g'
///
void MultiplyPointwiseRR(const grid_t *restrict_ f, const grid_t *restrict_ g, grid_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < N_GRID*N_GRID; i++)
	{
		ret->data[i] = f->data[i] * g->data[i];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'g'
///
void MultiplyPointwiseRC(const grid_t *restrict_ f, const cgrid_t *restrict_ g, cgrid_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < N_GRID*N_GRID; i++)
	{
		ret->data[i][0] = f->data[i] * g->data[i][0];
		ret->data[i][1] = f->data[i] * g->data[i][1];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'g'
///
void MultiplyPointwiseCC(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, cgrid_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < N_GRID*N_GRID; i++)
	{
		// complex multiplication
		ret->data[i][0] = f->data[i][0] * g->data[i][0] - f->data[i][1] * g->data[i][1];
		ret->data[i][1] = f->data[i][0] * g->data[i][1] + f->data[i][1] * g->data[i][0];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'g'
///
void MultiplyPointwiseRR2(const grid2_t *restrict_ f, const grid2_t *restrict_ g, grid2_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < 4*N_GRID*N_GRID; i++)
	{
		ret->data[i] = f->data[i] * g->data[i];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'g'
///
void MultiplyPointwiseCC2(const cgrid2_t *restrict_ f, const cgrid2_t *restrict_ g, cgrid2_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < 4*N_GRID*N_GRID; i++)
	{
		// complex multiplication
		ret->data[i][0] = f->data[i][0] * g->data[i][0] - f->data[i][1] * g->data[i][1];
		ret->data[i][1] = f->data[i][0] * g->data[i][1] + f->data[i][1] * g->data[i][0];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'ret', such that ret[i] *= f[i]
///
void MultiplyPointwiseAssignRR(const grid_t *restrict_ f, grid_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < N_GRID*N_GRID; i++)
	{
		ret->data[i] *= f->data[i];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'ret', such that ret[i] *= f[i]
///
void MultiplyPointwiseAssignRC(const grid_t *restrict_ f, cgrid_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < N_GRID*N_GRID; i++)
	{
		ret->data[i][0] *= f->data[i];
		ret->data[i][1] *= f->data[i];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'ret', such that ret[i] *= f[i]
///
void MultiplyPointwiseAssignRC2(const grid2_t *restrict_ f, cgrid2_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < 4*N_GRID*N_GRID; i++)
	{
		ret->data[i][0] *= f->data[i];
		ret->data[i][1] *= f->data[i];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Pointwise multiplication of the entries in grids 'f' and 'ret', such that ret[i] *= f[i]
///
void MultiplyPointwiseAssignCC2(const cgrid2_t *restrict_ f, cgrid2_t *restrict_ ret)
{
	unsigned int i;
	for (i = 0; i < 4*N_GRID*N_GRID; i++)
	{
		// complex multiplication
		// store current value of 'ret' since real part will be overwritten in the next line
		double real = ret->data[i][0];
		double imag = ret->data[i][1];
		ret->data[i][0] = f->data[i][0] * real - f->data[i][1] * imag;
		ret->data[i][1] = f->data[i][0] * imag + f->data[i][1] * real;
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Zero-pad grid to double the size in each direction
///
void GridZeroPad(const cgrid_t *restrict_ f, cgrid2_t *restrict_ ret)
{
	unsigned int l;

	memset(ret->data, 0, sizeof(ret->data));

	for (l = 0; l < N_GRID / 2; l++)
	{
		memcpy(&ret->data[             2*N_GRID*l], &f->data[           N_GRID*l], (N_GRID/2)*sizeof(fftw_complex));
		memcpy(&ret->data[3*N_GRID/2 + 2*N_GRID*l], &f->data[N_GRID/2 + N_GRID*l], (N_GRID/2)*sizeof(fftw_complex));
	}

	for (l = N_GRID / 2; l < N_GRID; l++)
	{
		memcpy(&ret->data[             2*N_GRID*(N_GRID + l)], &f->data[           N_GRID*l], (N_GRID/2)*sizeof(fftw_complex));
		memcpy(&ret->data[3*N_GRID/2 + 2*N_GRID*(N_GRID + l)], &f->data[N_GRID/2 + N_GRID*l], (N_GRID/2)*sizeof(fftw_complex));
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Extract the inner quarter of the larger grid (reverses zero-padding)
///
void GridExtract(const cgrid2_t *restrict_ f, cgrid_t *restrict_ ret)
{
	unsigned int l;

	for (l = 0; l < N_GRID / 2; l++)
	{
		memcpy(&ret->data[           N_GRID*l], &f->data[             2*N_GRID*l], (N_GRID/2)*sizeof(fftw_complex));
		memcpy(&ret->data[N_GRID/2 + N_GRID*l], &f->data[3*N_GRID/2 + 2*N_GRID*l], (N_GRID/2)*sizeof(fftw_complex));
	}

	for (l = N_GRID / 2; l < N_GRID; l++)
	{
		memcpy(&ret->data[           N_GRID*l], &f->data[             2*N_GRID*(N_GRID + l)], (N_GRID/2)*sizeof(fftw_complex));
		memcpy(&ret->data[N_GRID/2 + N_GRID*l], &f->data[3*N_GRID/2 + 2*N_GRID*(N_GRID + l)], (N_GRID/2)*sizeof(fftw_complex));
	}
}
