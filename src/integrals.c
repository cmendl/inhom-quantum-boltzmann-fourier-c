/// \file integrals.c
/// \brief Evaluation of the I1 (and I4), I2 and I3 integrals given a quadrature formula representing the delta function or principal value
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

#include "integrals.h"
#include "util.h"
#include <memory.h>



//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data and create FFTW plans for I1 integral
///
void I1interm_Create(I1interm_t* interm)
{
	interm->f2  = fftw_malloc(sizeof(interm->f2[0]));
	interm->f2F = fftw_malloc(sizeof(interm->f2F[0]));

	interm->h2  = fftw_malloc(sizeof(interm->h2[0]));
	interm->h2F = fftw_malloc(sizeof(interm->h2F[0]));

	interm->tmp = fftw_malloc(sizeof(interm->tmp[0]));

	interm->tmp2  = fftw_malloc(sizeof(interm->tmp2[0]));
	interm->tmp2F = fftw_malloc(sizeof(interm->tmp2F[0]));

	interm->plan_f    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->f2->data,    interm->f2F->data,   FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_h    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->h2->data,    interm->h2F->data,   FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_forw = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->tmp2->data,  interm->tmp2F->data, FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_back = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->tmp2F->data, interm->tmp2->data,  FFTW_BACKWARD, FFTW_ESTIMATE);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data and destroy FFTW plans used for I1 integral
///
void I1interm_Delete(I1interm_t* interm)
{
	fftw_destroy_plan(interm->plan_back);
	fftw_destroy_plan(interm->plan_forw);
	fftw_destroy_plan(interm->plan_h);
	fftw_destroy_plan(interm->plan_f);

	fftw_free(interm->tmp2F);
	fftw_free(interm->tmp2);
	interm->tmp2F = NULL;
	interm->tmp2  = NULL;

	fftw_free(interm->tmp);
	interm->tmp = NULL;

	fftw_free(interm->h2F);
	fftw_free(interm->h2);
	interm->h2F = NULL;
	interm->h2  = NULL;

	fftw_free(interm->f2F);
	fftw_free(interm->f2);
	interm->f2F = NULL;
	interm->f2  = NULL;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate I1 integral
///
void I1integral(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, const cgrid_t *restrict_ h, const quadI1_t *restrict_ quad, I1interm_t *restrict_ interm, cgrid_t *restrict_ I1)
{
	unsigned int i, j;
	const double scale = square(1.0 / (4 * N_GRID * N_GRID));

	// precompute FFT of extended f
	GridZeroPad(f, interm->f2);
	fftw_execute(interm->plan_f);

	// precompute FFT of extended h
	GridZeroPad(h, interm->h2);
	fftw_execute(interm->plan_h);

	memset(&I1->data, 0, sizeof(I1->data));

	for (j = 0; j < quad->num; j++)
	{
		MultiplyPointwiseRC(&quad->psiR2[j], g, interm->tmp);
		GridZeroPad(interm->tmp, interm->tmp2);
		fftw_execute(interm->plan_forw);		// FFT of zero-padded g_psi

		// convolution of f with g_psi
		MultiplyPointwiseAssignCC2(interm->f2F, interm->tmp2F);
		fftw_execute(interm->plan_back);

		// pointwise multiplication with psiR1
		MultiplyPointwiseAssignRC2(&quad->psiR1[j], interm->tmp2);

		// convolution with h
		fftw_execute(interm->plan_forw);
		MultiplyPointwiseAssignCC2(interm->h2F, interm->tmp2F);
		fftw_execute(interm->plan_back);

		// extract [-N/2:N/2-1] part
		GridExtract(interm->tmp2, interm->tmp);

		// could use LAPACK here...?
		for (i = 0; i < N_GRID*N_GRID; i++)
		{
			I1->data[i][0] += scale * interm->tmp->data[i][0];
			I1->data[i][1] += scale * interm->tmp->data[i][1];
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data and create FFTW plans for I2 integral
///
void I2interm_Create(I2interm_t* interm)
{
	interm->f_psi = fftw_malloc(sizeof(interm->f_psi[0]));
	interm->g_psi = fftw_malloc(sizeof(interm->g_psi[0]));

	interm->f_psi2 = fftw_malloc(sizeof(interm->f_psi2[0]));
	interm->g_psi2 = fftw_malloc(sizeof(interm->g_psi2[0]));

	interm->f_psi2F = fftw_malloc(sizeof(interm->f_psi2F[0]));
	interm->g_psi2F = fftw_malloc(sizeof(interm->g_psi2F[0]));

	interm->h2  = fftw_malloc(sizeof(interm->h2[0]));
	interm->h2F = fftw_malloc(sizeof(interm->h2F[0]));

	interm->tmp = fftw_malloc(sizeof(interm->tmp[0]));

	interm->tmp2  = fftw_malloc(sizeof(interm->tmp2[0]));
	interm->tmp2F = fftw_malloc(sizeof(interm->tmp2F[0]));

	interm->plan_f    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->f_psi2->data, interm->f_psi2F->data, FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_g    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->g_psi2->data, interm->g_psi2F->data, FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_h    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->h2->data,     interm->h2F->data,     FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_back = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->tmp2F->data,  interm->tmp2->data,    FFTW_BACKWARD, FFTW_ESTIMATE);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data and destroy FFTW plans used for I2 integral
///
void I2interm_Delete(I2interm_t* interm)
{
	fftw_destroy_plan(interm->plan_back);
	fftw_destroy_plan(interm->plan_h);
	fftw_destroy_plan(interm->plan_g);
	fftw_destroy_plan(interm->plan_f);

	fftw_free(interm->tmp2F);
	fftw_free(interm->tmp2);
	interm->tmp2F = NULL;
	interm->tmp2  = NULL;

	fftw_free(interm->tmp);
	interm->tmp = NULL;

	fftw_free(interm->h2F);
	fftw_free(interm->h2);
	interm->h2F = NULL;
	interm->h2  = NULL;

	fftw_free(interm->g_psi2F);
	fftw_free(interm->f_psi2F);
	interm->g_psi2F = NULL;
	interm->f_psi2F = NULL;

	fftw_free(interm->g_psi2);
	fftw_free(interm->f_psi2);
	interm->g_psi2 = NULL;
	interm->f_psi2 = NULL;

	fftw_free(interm->g_psi);
	fftw_free(interm->f_psi);
	interm->g_psi = NULL;
	interm->f_psi = NULL;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate I2 integral
///
void I2integral(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, const cgrid_t *restrict_ h, const quadI2_t *restrict_ quad, I2interm_t *restrict_ interm, cgrid_t *restrict_ I2)
{
	unsigned int i, j;
	const double scale = 1.0 / (4 * N_GRID * N_GRID);

	// precompute FFT of extended h
	GridZeroPad(h, interm->h2);
	fftw_execute(interm->plan_h);

	memset(&I2->data, 0, sizeof(I2->data));

	for (j = 0; j < quad->num; j++)
	{
		MultiplyPointwiseRC(&quad->psiR1[j], f, interm->f_psi);
		MultiplyPointwiseRC(&quad->psiR2[j], g, interm->g_psi);

		// zero-pad
		GridZeroPad(interm->f_psi, interm->f_psi2);
		GridZeroPad(interm->g_psi, interm->g_psi2);

		// FFT of f_psi2 and g_psi2
		fftw_execute(interm->plan_f);
		fftw_execute(interm->plan_g);

		// convolution of f_psi2, g_psi2 and h2F
		// accumulate products of f_psi2F, g_psi2F and h2F in tmp2F
		MultiplyPointwiseCC2(interm->f_psi2F, interm->g_psi2F, interm->tmp2F);
		MultiplyPointwiseAssignCC2(interm->h2F, interm->tmp2F);
		fftw_execute(interm->plan_back);

		// extract [-N/2:N/2-1] part
		GridExtract(interm->tmp2, interm->tmp);

		// could use LAPACK here...?
		for (i = 0; i < N_GRID*N_GRID; i++)
		{
			I2->data[i][0] += scale * interm->tmp->data[i][0];
			I2->data[i][1] += scale * interm->tmp->data[i][1];
		}
	}
}



//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data and create FFTW plans for I3 integral
///
void I3interm_Create(I3interm_t* interm)
{
	interm->f_psi = fftw_malloc(sizeof(interm->f_psi[0]));
	interm->h_psi = fftw_malloc(sizeof(interm->h_psi[0]));

	interm->f_psi2 = fftw_malloc(sizeof(interm->f_psi2[0]));
	interm->h_psi2 = fftw_malloc(sizeof(interm->h_psi2[0]));

	interm->f_psi2F = fftw_malloc(sizeof(interm->f_psi2F[0]));
	interm->h_psi2F = fftw_malloc(sizeof(interm->h_psi2F[0]));

	interm->g2  = fftw_malloc(sizeof(interm->g2[0]));
	interm->g2F = fftw_malloc(sizeof(interm->g2F[0]));

	interm->tmp = fftw_malloc(sizeof(interm->tmp[0]));

	interm->tmp2  = fftw_malloc(sizeof(interm->tmp2[0]));
	interm->tmp2F = fftw_malloc(sizeof(interm->tmp2F[0]));

	interm->plan_f    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->f_psi2->data, interm->f_psi2F->data, FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_h    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->h_psi2->data, interm->h_psi2F->data, FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_g    = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->g2->data,     interm->g2F->data,     FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_forw = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->tmp2->data,   interm->tmp2F->data,   FFTW_FORWARD,  FFTW_ESTIMATE);
	interm->plan_back = fftw_plan_dft_2d(2*N_GRID, 2*N_GRID, interm->tmp2F->data,  interm->tmp2->data,    FFTW_BACKWARD, FFTW_ESTIMATE);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data and destroy FFTW plans used for I3 integral
///
void I3interm_Delete(I3interm_t* interm)
{
	fftw_destroy_plan(interm->plan_back);
	fftw_destroy_plan(interm->plan_forw);
	fftw_destroy_plan(interm->plan_g);
	fftw_destroy_plan(interm->plan_h);
	fftw_destroy_plan(interm->plan_f);

	fftw_free(interm->tmp2F);
	fftw_free(interm->tmp2);
	interm->tmp2F = NULL;
	interm->tmp2  = NULL;

	fftw_free(interm->tmp);
	interm->tmp = NULL;

	fftw_free(interm->g2F);
	fftw_free(interm->g2);
	interm->g2F = NULL;
	interm->g2  = NULL;

	fftw_free(interm->h_psi2F);
	fftw_free(interm->f_psi2F);
	interm->h_psi2F = NULL;
	interm->f_psi2F = NULL;

	fftw_free(interm->h_psi2);
	fftw_free(interm->f_psi2);
	interm->h_psi2 = NULL;
	interm->f_psi2 = NULL;

	fftw_free(interm->h_psi);
	fftw_free(interm->f_psi);
	interm->h_psi = NULL;
	interm->f_psi = NULL;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate I3 integral
///
void I3integral(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, const cgrid_t *restrict_ h, const quadI3_t *restrict_ quad, I3interm_t *restrict_ interm, cgrid_t *restrict_ I3)
{
	unsigned int i, j;
	const double scale = square(1.0 / (4 * N_GRID * N_GRID));

	// precompute FFT of extended g
	GridZeroPad(g, interm->g2);
	fftw_execute(interm->plan_g);

	memset(&I3->data, 0, sizeof(I3->data));

	for (j = 0; j < quad->num; j++)
	{
		// pointwise multiplication of two matrices
		MultiplyPointwiseCC(&quad->psiR1[j], f, interm->f_psi);
		MultiplyPointwiseCC(&quad->psiR1[j], h, interm->h_psi);

		// zero-padding
		GridZeroPad(interm->f_psi, interm->f_psi2);
		GridZeroPad(interm->h_psi, interm->h_psi2);

		// FFT of f_psi2 and h_psi2
		fftw_execute(interm->plan_f);
		fftw_execute(interm->plan_h);

		// convolution of g2 and h_psi2
		MultiplyPointwiseCC2(interm->g2F, interm->h_psi2F, interm->tmp2F);
		fftw_execute(interm->plan_back);

		// pointwise multiplication with psiR2[j]
		MultiplyPointwiseAssignRC2(&quad->psiR2[j], interm->tmp2);

		// convolution with f_psi
		fftw_execute(interm->plan_forw);
		MultiplyPointwiseAssignCC2(interm->f_psi2F, interm->tmp2F);
		fftw_execute(interm->plan_back);

		// extract [-N/2:N/2-1] part
		GridExtract(interm->tmp2, interm->tmp);

		// could use LAPACK here...?
		for (i = 0; i < N_GRID*N_GRID; i++)
		{
			I3->data[i][0] += scale * interm->tmp->data[i][0];
			I3->data[i][1] += scale * interm->tmp->data[i][1];
		}
	}
}


// 'I4integral' is the same as 'I1integral'
