/// \file quadrature.c
/// \brief Construction of quadrature formulas representing the Fourier-transformed delta function or principal value, as input for the I1, I2, I3 and I4 integrals
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

#include "quadrature.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <memory.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for I1 quadrature structure
///
void QuadI1_Allocate(const unsigned int num, quadI1_t *quad)
{
	// use 'fftw_malloc' for proper memory alignment
	quad->num = num;
	quad->psiR1 = fftw_malloc(num * sizeof(quad->psiR1[0]));
	quad->psiR2 = fftw_malloc(num * sizeof(quad->psiR2[0]));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete I1 quadrature structure (free memory)
///
void QuadI1_Delete(quadI1_t *quad)
{
	fftw_free(quad->psiR2);
	fftw_free(quad->psiR1);
	quad->num = 0;
}



//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for I2 quadrature structure
///
void QuadI2_Allocate(const unsigned int num, quadI2_t *quad)
{
	// use 'fftw_malloc' for proper memory alignment
	quad->num = num;
	quad->psiR1 = fftw_malloc(num * sizeof(quad->psiR1[0]));
	quad->psiR2 = fftw_malloc(num * sizeof(quad->psiR2[0]));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete I2 quadrature structure (free memory)
///
void QuadI2_Delete(quadI2_t *quad)
{
	fftw_free(quad->psiR2);
	fftw_free(quad->psiR1);
	quad->num = 0;
}



//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for I3 quadrature structure
///
void QuadI3_Allocate(const unsigned int num, quadI3_t *quad)
{
	// use 'fftw_malloc' for proper memory alignment
	quad->num = num;
	quad->psiR1 = fftw_malloc(num * sizeof(quad->psiR1[0]));
	quad->psiR2 = fftw_malloc(num * sizeof(quad->psiR2[0]));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete I3 quadrature structure (free memory)
///
void QuadI3_Delete(quadI3_t *quad)
{
	fftw_free(quad->psiR2);
	fftw_free(quad->psiR1);
	quad->num = 0;
}



//_______________________________________________________________________________________________________________________
///
/// \brief sin(x) * sinc(x)
///
static inline double sin_sinc(double x)
{
	if (fabs(x) < 1e-8)
	{
		return x;
	}
	else
	{
		return square(sin(x))/x;
	}
}

//_______________________________________________________________________________________________________________________
///
/// \brief sinc(x)
///
static inline double sinc(double x)
{
	if (fabs(x) < 1e-8)
	{
		return 1;
	}
	else
	{
		return sin(x)/x;
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Fourier transform of principal value of dot product between to vectors u and u' used in I1 integral,
/// G(xi,chi) = int_{B_R} du int_{B_R} du' P(1/(u cdot u')) exp(i pi xi cdot u/L) exp(i pi chi cdot u'/L)
///
void FourierI1(const unsigned int J, const double L, const double R, quadI1_t *quad)
{
	unsigned int i, j, k, l;

	double *theta1 = malloc(J * sizeof(double));
	double *theta2 = malloc(J * sizeof(double));
	const double deltaTh = M_PI / J;
	for (j = 0; j < J; j++)
	{
		theta1[j] = j*deltaTh;
		theta2[j] = theta1[j] + 0.5*deltaTh;
	}

	int kgrid[N_GRID];
	for (k = 0; k < N_GRID; k++) {
		kgrid[k] = k < N_GRID/2 ? k : k - N_GRID;
	}
	int kgrid2[2*N_GRID];
	for (k = 0; k < 2*N_GRID; k++) {
		kgrid2[k] = k < N_GRID ? k : k - 2*N_GRID;
	}

	// allocate memory
	QuadI1_Allocate(J*J, quad);

	const double s = R * M_PI/(2*L);

	for (i = 0; i < J; i++)
	{
		double cos_th1 = cos(theta1[i]);
		double sin_th1 = sin(theta1[i]);

		for (l = 0; l < 2*N_GRID; l++)
		{
			for (k = 0; k < 2*N_GRID; k++)
			{
				// follow convention of 'meshgrid' in Matlab, to agree with Matlab implementation for J odd
				quad->psiR1[i].data[k + 2*N_GRID*l] = R * sin_sinc(s * (kgrid2[l]*cos_th1 + kgrid2[k]*sin_th1));
			}
		}

		// duplicate
		for (j = 1; j < J; j++)		// loop starts at 1!
		{
			memcpy(&quad->psiR1[i + J*j], &quad->psiR1[i], sizeof(quad->psiR1[i]));
		}
	}

	for (j = 0; j < J; j++)
	{
		double cos_th2 = cos(theta2[j]);
		double sin_th2 = sin(theta2[j]);

		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// follow convention of 'meshgrid' in Matlab, to agree with Matlab implementation for J odd
				quad->psiR2[J*j].data[k + N_GRID*l] = R * sin_sinc(s * (kgrid[l]*cos_th2 + kgrid[k]*sin_th2));
			}
		}

		// duplicate
		for (i = 1; i < J; i++)		// loop starts at 1!
		{
			memcpy(&quad->psiR2[i + J*j], &quad->psiR2[J*j], sizeof(quad->psiR2[J*j]));
		}
	}

	for (j = 0; j < J; j++)
	{
		for (i = 0; i < J; i++)
		{
			// principal value singularity, minus sign due to i^2
			double weight = -square(2*M_PI/J)/cos(theta1[i] - theta2[j]);

			for (l = 0; l < N_GRID; l++)
			{
				for (k = 0; k < N_GRID; k++)
				{
					quad->psiR2[i + J*j].data[k + N_GRID*l] *= weight;
				}
			}
		}
	}

	free(theta2);
	free(theta1);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Fourier transform of delta function of dot product between to vectors used in I2 integral
///
void FourierI2(const unsigned int J, const double L, const double R, quadI2_t *quad)
{
	unsigned int j, k, l;

	double *theta = malloc(J * sizeof(double));
	const double deltaTh = M_PI / J;
	for (j = 0; j < J; j++)
	{
		theta[j] = j*deltaTh;
	}

	int kgrid[N_GRID];
	for (k = 0; k < N_GRID; k++) {
		kgrid[k] = k < N_GRID/2 ? k : k - N_GRID;
	}

	// allocate memory
	QuadI2_Allocate(J, quad);

	const double f1 = 2*R;
	const double f2 = 2*R*M_PI/J;	// factor pi/J due to uniform angular integration weight

	const double s = R * M_PI/L;

	for (j = 0; j < J; j++)
	{
		double cos_th = cos(theta[j]);
		double sin_th = sin(theta[j]);

		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// follow convention of 'meshgrid' in Matlab
				quad->psiR1[j].data[k + N_GRID*l] = f1 * sinc(s * ( kgrid[l]*cos_th + kgrid[k]*sin_th));
				quad->psiR2[j].data[k + N_GRID*l] = f2 * sinc(s * (-kgrid[l]*sin_th + kgrid[k]*cos_th));	// shift by pi/2
			}
		}
	}

	free(theta);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Compute the Legendre-Gauss nodes and weights on an interval [a,b] with truncation order M-1
///
void LegendreGaussQuad(const unsigned int M, const double a, const double b, double *x, double *w)
{
	// This script is for computing definite integrals using Legendre-Gauss
	// Quadrature. Computes the Legendre-Gauss nodes and weights on an interval
	// [a,b] with truncation order M-1
	//
	// Suppose you have a continuous function f(x) which is defined on [a,b]
	// which you can evaluate at any x in [a,b]. Simply evaluate it at all of
	// the values contained in the x vector to obtain a vector f. Then compute
	// the definite integral using sum(f.*w);
	//
	// Written by Greg von Winckel - 02 / 25 / 2004

	unsigned int i, k;

	const double invN = 1.0 / (M-1);

	double *xu = malloc(M * sizeof(double));
	for (i = 0; i < M; i++)
	{
		xu[i] = -1.0 + 2*i*invN;
	}

	// initial guess
	double *y = malloc(M * sizeof(double));
	for (i = 0; i < M; i++)
	{
		y[i] = cos((2*i+1)*M_PI/(2*M)) + 0.27/M * sin(M_PI*xu[i]*(M-1)/(M+1));
	}

	// Legendre-Gauss Vandermonde matrix (LGVM) with dimension (N+1) x (N+2)
	double *L = calloc(M*(M+1), sizeof(double));

	// derivative of LGVM
	double *Lp = calloc(M, sizeof(double));

	// compute the zeros of the M-th Legendre polynomial
	// using the recursion relation and the Newton-Raphson method

	double *y0 = malloc(M * sizeof(double));
	for (i = 0; i < M; i++) {
		y0[i] = 2;
	}

	// iterate until new points are uniformly within epsilon of old points
	double diff = 0;
	for (i = 0; i < M; i++) {
		diff = maxf(diff, fabs(y[i] - y0[i]));
	}
	while (diff > DBL_EPSILON)
	{
		for (i = 0; i < M; i++)
		{
			L[i] = 1;
			L[i + M] = y[i];
		}

		for (k = 1; k < M; k++)
		{
			for (i = 0; i < M; i++) {
				L[i +  (k + 1)*M] = ((2*k + 1)*y[i] * L[i + k*M] - k*L[i + (k - 1)*M]) / (k + 1);
			}
		}

		for (i = 0; i < M; i++) {
			Lp[i] = (M+1)*(L[i + (M-1)*M] - y[i]*L[i + M*M]) / (1 - square(y[i]));
		}

		memcpy(y0, y, M*sizeof(double));

		diff = 0;
		for (i = 0; i < M; i++)
		{
			double t = L[i + M*M] / Lp[i];
			y[i] -= t;

			diff = maxf(diff, fabs(t));
		}
	}

	for (i = 0; i < M; i++)
	{
		// linear map from [-1,1] to [a,b]
		x[i] = (a*(1 - y[i]) + b*(1 + y[i])) / 2;

		// compute the weights
		w[i] = (b - a) / ((1 - square(y[i]))*square(Lp[i]))*square((double)(M+1) / M);
	}

	// clean up
	free(L);
	free(Lp);
	free(y0);
	free(y);
	free(xu);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Fourier transform of delta function of dot product between to vectors used in I3 integral
///
void FourierI3(const unsigned int J, const double L, const double R, const unsigned int M, quadI3_t *quad)
{
	unsigned int j, m, k, l;

	// uniform discretization of the interval [0,pi]
	double *theta = malloc(J * sizeof(double));
	const double deltaTh = M_PI / J;
	for (j = 0; j < J; j++)
	{
		theta[j] = j*deltaTh;
	}

	int kgrid[N_GRID];
	for (k = 0; k < N_GRID; k++) {
		kgrid[k] = k < N_GRID/2 ? k : k - N_GRID;
	}
	int kgrid2[2*N_GRID];
	for (k = 0; k < 2*N_GRID; k++) {
		kgrid2[k] = k < N_GRID ? k : k - 2*N_GRID;
	}

	// allocate memory
	QuadI3_Allocate(J * M, quad);

	// Gauss-Legendre quadrature
	double *rho = malloc(M * sizeof(double));
	double *w   = malloc(M * sizeof(double));
	LegendreGaussQuad(M, -R, R, rho, w);

	const double s2 = M_PI * R/L;

	for (j = 0; j < J; j++)
	{
		double cos_th = cos(theta[j]);
		double sin_th = sin(theta[j]);

		for (m = 0; m < M; m++)
		{
			cgrid_t *curgrid1 = &quad->psiR1[m + j*M];
			grid2_t *curgrid2 = &quad->psiR2[m + j*M];

			const double s1 = rho[m] * M_PI/L;

			// include integration weights for radial and angular integral
			const double f2 = 2*R * w[m] * M_PI/J;

			for (l = 0; l < N_GRID; l++)
			{
				for (k = 0; k < N_GRID; k++)
				{
					double t = s1 * (kgrid[l]*cos_th + kgrid[k]*sin_th);	// follow convention of 'meshgrid' in Matlab, to agree with Matlab implementation for J odd
					curgrid1->data[k + N_GRID*l][0] = cos(t);
					curgrid1->data[k + N_GRID*l][1] = sin(t);
				}
			}

			for (l = 0; l < 2*N_GRID; l++)
			{
				for (k = 0; k < 2*N_GRID; k++)
				{
					// shift by pi/2
					// follow convention of 'meshgrid' in Matlab, to agree with Matlab implementation for J odd
					curgrid2->data[k + 2*N_GRID*l] = f2 * sinc(s2 * (-kgrid2[l]*sin_th + kgrid2[k]*cos_th));
				}
			}
		}
	}

	free(w);
	free(rho);
	free(theta);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Fourier transform of delta function of dot product between to vectors used in I4 integral
///
void FourierI4(const unsigned int J, const double L, const double R, quadI1_t *quad)
{
	unsigned int j, k, l;

	// uniform discretization of the interval [0,pi]
	double *theta = malloc(J * sizeof(double));
	const double deltaTh = M_PI / J;
	for (j = 0; j < J; j++)
	{
		theta[j] = j*deltaTh;
	}

	int kgrid[N_GRID];
	for (k = 0; k < N_GRID; k++) {
		kgrid[k] = k < N_GRID/2 ? k : k - N_GRID;
	}
	int kgrid2[2*N_GRID];
	for (k = 0; k < 2*N_GRID; k++) {
		kgrid2[k] = k < N_GRID ? k : k - 2*N_GRID;
	}

	// allocate memory
	QuadI1_Allocate(J, quad);

	const double f1 = 2*R;
	const double f2 = 2*R*M_PI/J;	// factor pi/J due to uniform angular integration weight

	const double s = M_PI * R/L;

	for (j = 0; j < J; j++)
	{
		double cos_th = cos(theta[j]);
		double sin_th = sin(theta[j]);

		for (l = 0; l < 2*N_GRID; l++)
		{
			for (k = 0; k < 2*N_GRID; k++)
			{
				// follow convention of 'meshgrid' in Matlab, to agree with Matlab implementation for J odd
				quad->psiR1[j].data[k + 2*N_GRID*l] = f1 * sinc(s * (kgrid2[l]*cos_th + kgrid2[k]*sin_th));
			}
		}

		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// follow convention of 'meshgrid' in Matlab, to agree with Matlab implementation for J odd
				quad->psiR2[j].data[k + N_GRID*l] = f2 * sinc(s * (-kgrid[l]*sin_th + kgrid[k]*cos_th));	// shift by pi/2
			}
		}
	}

	free(theta);
}
