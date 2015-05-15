/// \file collision.c
/// \brief Calculate the collision operator using Fourier transformations
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

#include "collision.h"
#include <memory.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data for Cd collision operator
///
void CdInterm_Create(CdInterm_t *interm)
{
	I2interm_Create(&interm->intermI2);
	I3interm_Create(&interm->intermI3);
	I1interm_Create(&interm->intermI4);

	interm->Wt  = fftw_malloc(sizeof(wignerF_t));
	interm->tmp = fftw_malloc(sizeof(interm->tmp[0]));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data for Cd collision operator
///
void CdInterm_Delete(CdInterm_t *interm)
{
	fftw_free(interm->tmp);
	fftw_free(interm->Wt);
	interm->tmp = NULL;
	interm->Wt  = NULL;

	I1interm_Delete(&interm->intermI4);
	I3interm_Delete(&interm->intermI3);
	I2interm_Delete(&interm->intermI2);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Collision integral of dissipative collision operator Cd
///
/// W is a 1 x 4 cell of N x N matrices, representing the Fourier-transformed
/// Wigner matrices in the Pauli-spin basis (including identity as first component).
/// The following notations are used:
///
/// W_i = w_{i,tr} id + w_{i,x} sigma_x + w_{i,y} sigma_y + w_{i,z} sigma_z
///     = w_i cdot sigma
///
/// eta = diag(1,-1,-1,-1)
///
void CdInt(const wignerF_t *restrict_ W, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4, CdInterm_t *restrict_ interm, wignerF_t *restrict_ Cd)
{
	unsigned int i, j, k;
	const double  eta[4] = { 1, -1, -1, -1 };
	const double eta2[4] = { 2, -2, -2, -2 };

	memset(Cd, 0, sizeof(Cd[0]));

	// tilde{W} = 1 - W in Fourier space
	for (j = 0; j < 4; j++)
	{
		for (k = 0; k < N_GRID*N_GRID; k++)
		{
			interm->Wt->comp[j].data[k][0] = -W->comp[j].data[k][0];
			interm->Wt->comp[j].data[k][1] = -W->comp[j].data[k][1];
		}
	}
	// Fourier transform of constant N x N matrix with all entries 1 is zero except for entry corresponding to frequency zero
	interm->Wt->comp[0].data[0][0]++;

	// 2 <w_3, w_4>_{\eta} (id - W_1)
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 4; i++)
		{
			I2integral(&W->comp[i], &W->comp[i], &interm->Wt->comp[j], quadI2, &interm->intermI2, interm->tmp);
			for (k = 0; k < N_GRID*N_GRID; k++)
			{
				Cd->comp[j].data[k][0] += eta2[i] * interm->tmp->data[k][0];
				Cd->comp[j].data[k][1] += eta2[i] * interm->tmp->data[k][1];
			}
		}
	}


	// - 2 <w_3, w_4>_{\eta} (\eta w_2)\cdot\sigma
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 4; i++)
		{
			I3integral(&W->comp[i], &W->comp[i], &W->comp[j], quadI3, &interm->intermI3, interm->tmp);
			for (k = 0; k < N_GRID*N_GRID; k++)
			{
				Cd->comp[j].data[k][0] -= eta[i] * eta2[j] * interm->tmp->data[k][0];
				Cd->comp[j].data[k][1] -= eta[i] * eta2[j] * interm->tmp->data[k][1];
			}
		}
	}


	// (w_{3,tr} + w_{3,tr} - 1) anticomm( (\eta w_2)\cdot\sigma, W_1 ) =
	// 2 (2*w_{3,tr} - 1) ( (w_{1,tr}(eta w_2) + w_{2,tr} w_1)\cdot\sigma - <w_1, w_2>_{id} id)

	// 2*W_0 - 1 in Fourier representation
	for (k = 0; k < N_GRID*N_GRID; k++)
	{
		interm->Wt->comp[0].data[k][0] = 2 * W->comp[0].data[k][0];
		interm->Wt->comp[0].data[k][1] = 2 * W->comp[0].data[k][1];
	}
	// Fourier transform of constant N x N matrix with all entries 1 is zero except for entry corresponding to frequency zero
	interm->Wt->comp[0].data[0][0]--;

	// sub-expression 2 (2*w_{3,tr} - 1) (w_{1,tr}(eta w_2))\cdot\sigma
	for (j = 0; j < 4; j++)
	{
		I1integral(&interm->Wt->comp[0], &W->comp[j], &W->comp[0], quadI4, &interm->intermI4, interm->tmp);
		for (k = 0; k < N_GRID*N_GRID; k++)
		{
			Cd->comp[j].data[k][0] += eta2[j] * interm->tmp->data[k][0];
			Cd->comp[j].data[k][1] += eta2[j] * interm->tmp->data[k][1];
		}
	}

	// sub-expression 2 (2*w_{3,tr} - 1) (w_{2,tr} w_1)\cdot\sigma
	for (j = 0; j < 4; j++)
	{
		I1integral(&interm->Wt->comp[0], &W->comp[0], &W->comp[j], quadI4, &interm->intermI4, interm->tmp);
		for (k = 0; k < N_GRID*N_GRID; k++)
		{
			Cd->comp[j].data[k][0] += 2 * interm->tmp->data[k][0];
			Cd->comp[j].data[k][1] += 2 * interm->tmp->data[k][1];
		}
	}

	// sub-expression - 2 (2*w_{3,tr} - 1) <w_1, w_2>_{id} id
	for (j = 0; j < 4; j++)
	{
		I1integral(&interm->Wt->comp[0], &W->comp[j], &W->comp[j], quadI4, &interm->intermI4, interm->tmp);
		for (k = 0; k < N_GRID*N_GRID; k++)
		{
			Cd->comp[0].data[k][0] -= 2 * interm->tmp->data[k][0];
			Cd->comp[0].data[k][1] -= 2 * interm->tmp->data[k][1];
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data for Cc collision operator
///
void CcInterm_Create(CcInterm_t *interm)
{
	I1interm_Create(&interm->intermI1);
	interm->Wtr2m1 = fftw_malloc(sizeof(interm->Wtr2m1[0]));
	interm->I1     = fftw_malloc(sizeof(interm->I1[0]));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data for Cc collision operator
///
void CcInterm_Delete(CcInterm_t *interm)
{
	fftw_free(interm->I1);
	fftw_free(interm->Wtr2m1);
	interm->I1     = NULL;
	interm->Wtr2m1 = NULL;

	I1interm_Delete(&interm->intermI1);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Collision integral of conservative collision operator Cc
///
/// W is a 1 x 4 cell of N x N matrices, representing the Fourier-transformed
/// Wigner matrices in the Pauli-spin basis (including identity as first component).
/// The following notations are used:
///
/// W_i = w_{i,tr} id + w_{i,x} sigma_x + w_{i,y} sigma_y + w_{i,z} sigma_z
///     = w_i cdot sigma
///
void CcInt(const wignerF_t *restrict_ W, const quadI1_t *restrict_ quadI1, CcInterm_t *restrict_ interm, wignerF_t *restrict_ Cc)
{
	unsigned int k;

	memset(Cc, 0, sizeof(Cc[0]));

	// trace terms of Cc are zero due to commutator

	// (w_{3,tr} + w_{3,tr} - 1) i [ W_2, W_1 ] =
	// - 2 (2*w_{3,tr} - 1) ( vec{w_2} cross vec{w_1} ) cdot sigma

	// -2 * (2*w_tr - 1)
	for (k = 0; k < N_GRID*N_GRID; k++)
	{
		interm->Wtr2m1->data[k][0] = -4*W->comp[0].data[k][0];
		interm->Wtr2m1->data[k][1] = -4*W->comp[0].data[k][1];
	}
	// Fourier transform of constant N x N matrix with all entries 1 is zero except for entry corresponding to frequency zero
	interm->Wtr2m1->data[0][0] += 2;

	// cross product

	// y, z
	I1integral(interm->Wtr2m1, &W->comp[2], &W->comp[3], quadI1, &interm->intermI1, &Cc->comp[1]);
	I1integral(interm->Wtr2m1, &W->comp[3], &W->comp[2], quadI1, &interm->intermI1, interm->I1);
	for (k = 0; k < N_GRID*N_GRID; k++) {
		Cc->comp[1].data[k][0] -= interm->I1->data[k][0];
		Cc->comp[1].data[k][1] -= interm->I1->data[k][1];
	}

	// z, x
	I1integral(interm->Wtr2m1, &W->comp[3], &W->comp[1], quadI1, &interm->intermI1, &Cc->comp[2]);
	I1integral(interm->Wtr2m1, &W->comp[1], &W->comp[3], quadI1, &interm->intermI1, interm->I1);
	for (k = 0; k < N_GRID*N_GRID; k++) {
		Cc->comp[2].data[k][0] -= interm->I1->data[k][0];
		Cc->comp[2].data[k][1] -= interm->I1->data[k][1];
	}

	// x, y
	I1integral(interm->Wtr2m1, &W->comp[1], &W->comp[2], quadI1, &interm->intermI1, &Cc->comp[3]);
	I1integral(interm->Wtr2m1, &W->comp[2], &W->comp[1], quadI1, &interm->intermI1, interm->I1);
	for (k = 0; k < N_GRID*N_GRID; k++) {
		Cc->comp[3].data[k][0] -= interm->I1->data[k][0];
		Cc->comp[3].data[k][1] -= interm->I1->data[k][1];
	}
}
