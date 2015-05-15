/// \file finite_volume.c
/// \brief Finite volume simulation using slope limiters with various boundary conditions
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

#include "finite_volume.h"
#include "util.h"
#include <stdlib.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief min-mod function, see Eq. (16.52) in Randall LeVeque. Numerical Methods for Conservation Laws (1992)
///
static inline double minmod(const double a, const double b)
{
	return 0.5 * (copysign(1, a) + copysign(1, b)) * minf(fabs(a), fabs(b));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Slope limiter step in one dimension with periodic boundary conditions
///
/// Numerically solve u_t + A u_x = 0 with periodic boundary conditions
/// using slope limiter method.
///
/// Un is a 1 x N matrix, with N the number of finite volumes
///
/// Reference: Randall LeVeque. Numerical Methods for Conservation Laws (1992)
///
void SlopeLimiterStepPeriodic(const double h, const double dt, const double A, const unsigned int N, const double *restrict_ Un, double *restrict_ Un1)
{
	unsigned int j;

	const double Aabs = fabs(A);

	// nu_p = dt*lambda_p/h
	const double nu = dt*A/h;

	// for slope-limiter term
	const double S    = 0.5 * A    * (copysign(1, nu) - nu);
	const double Sabs = 0.5 * Aabs * (copysign(1, nu) - nu);

	// Eq. (16.34)
	// alpha_j = R^{-1}*(U^n_{j+1} - U^n_j)
	double *alpha = malloc(N * sizeof(double));
	for (j = 0; j < N-1; j++) {
		alpha[j] = Un[j+1] - Un[j];
	}
	alpha[N-1] = Un[0] - Un[N-1];	// periodic

	// h * beta_j, Eq. (16.56)
	double *h_beta = malloc(N * sizeof(double));
	for (j = 1; j < N; j++) {
		h_beta[j] = minmod(alpha[j], alpha[j-1]);
	}
	h_beta[0] = minmod(alpha[0], alpha[N-1]);

	// first term in square brackets of Eq. (16.57) is basically Godunov step in Eq. (13.17)
	const double scale = dt/(2*h);
	for (j = 1; j < N-1; j++)
	{
		Un1[j] = Un[j] + scale * (
			- A*(    Un[j+1] -     Un[j-1]) + Aabs*(    Un[j+1] - 2*    Un[j] +     Un[j-1])
			- S*(h_beta[j+1] - h_beta[j-1]) + Sabs*(h_beta[j+1] - 2*h_beta[j] + h_beta[j-1]));
	}
	// j == 0
	{
		Un1[0] = Un[0] + scale * (
			- A*(    Un[1] -     Un[N-1]) + Aabs*(    Un[1] - 2*    Un[0] +     Un[N-1])
			- S*(h_beta[1] - h_beta[N-1]) + Sabs*(h_beta[1] - 2*h_beta[0] + h_beta[N-1]));
	}
	// j == N-1
	{
		Un1[N-1] = Un[N-1] + scale * (
			- A*(    Un[0] -     Un[N-2]) + Aabs*(    Un[0] - 2*    Un[N-1] +     Un[N-2])
			- S*(h_beta[0] - h_beta[N-2]) + Sabs*(h_beta[0] - 2*h_beta[N-1] + h_beta[N-2]));
	}

	// clean up
	free(h_beta);
	free(alpha);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Slope limiter step in one dimension with Dirichlet boundary conditions
///
/// Numerically solve u_t + A u_x = 0 with Dirichlet boundary conditions
/// using slope limiter method.
///
/// Un is a 1 x N matrix, with N the number of finite volumes
///
/// Reference: Randall LeVeque. Numerical Methods for Conservation Laws (1992)
///
void SlopeLimiterStepDirichlet(const double h, const double dt, const double A, const unsigned int N, const double *restrict_ Un, double *restrict_ Un1)
{
	unsigned int j;

	const double Aabs = fabs(A);

	// nu_p = dt*lambda_p/h
	const double nu = dt*A/h;

	// for slope-limiter term
	const double S    = 0.5 * A    * (copysign(1, nu) - nu);
	const double Sabs = 0.5 * Aabs * (copysign(1, nu) - nu);

	// Eq. (16.34)
	// alpha_j = R^{-1}*(U^n_{j+1} - U^n_j)
	double *alpha = malloc(N * sizeof(double));
	for (j = 0; j < N-1; j++) {
		alpha[j] = Un[j+1] - Un[j];
	}
	alpha[N-1] = 0;		// Dirichlet

	// h * beta_j, Eq. (16.56)
	double *h_beta = malloc(N * sizeof(double));
	for (j = 1; j < N-1; j++) {
		h_beta[j] = minmod(alpha[j], alpha[j-1]);
	}
	// minmod(a, 0) == 0, minmod(0, b) == 0
	h_beta[0]   = 0;
	h_beta[N-1] = 0;

	// first term in square brackets of Eq. (16.57) is basically Godunov step in Eq. (13.17)
	const double scale = dt/(2*h);
	for (j = 1; j < N-1; j++)
	{
		Un1[j] = Un[j] + scale * (
			- A*(    Un[j+1] -     Un[j-1]) + Aabs*(    Un[j+1] - 2*    Un[j] +     Un[j-1])
			- S*(h_beta[j+1] - h_beta[j-1]) + Sabs*(h_beta[j+1] - 2*h_beta[j] + h_beta[j-1]));
	}
	// Dirichlet boundary conditions
	if (A >= 0)
	{
		// substitute incoming part by given Dirichlet boundary states
		Un1[0] = Un[0];
		// outgoing part: pad copies of boundary states, i.e., replace index N by N-1
		Un1[N-1] = Un[N-1] + 2*scale * (
			- A*(    Un[N-1] -     Un[N-2])
			- S*(h_beta[N-1] - h_beta[N-2]));
	}
	else
	{
		// outgoing part: pad copies of boundary states, i.e., replace index -1 by 0
		Un1[0] = Un[0] + 2*scale * (
			- A*(    Un[1] -     Un[0])
			- S*(h_beta[1] - h_beta[0]));
		// substitute incoming part by given Dirichlet boundary states
		Un1[N-1] = Un[N-1];
	}

	// clean up
	free(h_beta);
	free(alpha);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate flux to the right
///
double CalculateRightwardFlux(const unsigned int m, const double *restrict_ A, const double *restrict_ U)
{
	unsigned int i;

	double flux = 0;
	for (i = 0; i < m/2; i++)
	{
		assert(A[i] >= 0);		// first half must only contain non-negative velocities
		flux += A[i] * U[i];
	}

	return flux;
}

//_______________________________________________________________________________________________________________________
///
/// \brief Calculate flux to the left
///
double CalculateLeftwardFlux(const unsigned int m, const double *restrict_ A, const double *restrict_ U)
{
	unsigned int i;

	double flux = 0;
	for (i = m/2; i < m; i++)
	{
		assert(A[i] <= 0);		// second half must only contain non-positive velocities
		flux += fabs(A[i]) * U[i];
	}

	return flux;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Slope limiter step in m dimensions with Maxwell boundary conditions
///
/// Numerically solve u_t + A u_x = 0 with Maxwell boundary conditions
/// using slope limiter method.
///
/// \param h spatial mesh width
/// \param dt time step
/// \param m number of entries at each spatial location
/// \param A transport velocities, vector of length m, must be arranged in decreasing order such that first half A[0,...,m/2-1] is non-negative and A[i] = -A[m-i-1]
/// \param N number of finite volumes
/// \param lambda accommodation coefficient for Maxwell reflection operator
/// \param fluxL predetermined total outgoing flux at left  boundary
/// \param fluxR predetermined total outgoing flux at right boundary
/// \param UmaxwL normalized left  boundary state for Maxwell boundary condition
/// \param UmaxwR normalized right boundary state for Maxwell boundary condition
/// \param Un input m x N matrix
/// \param Un1 output m x N matrix, assuming that memory has been allocated
///
/// Reference: Randall LeVeque. Numerical Methods for Conservation Laws (1992)
///
void SlopeLimiterStepMaxwell(const double h, const double dt, const unsigned int m, const double *restrict_ A, const unsigned int N,
	const double lambda, const double fluxL, const double fluxR, const double *restrict_ UmaxwL, const double *restrict_ UmaxwR, const double *restrict_ Un, double *restrict_ Un1)
{
	unsigned int i;
	int j;

	// Maxwell boundary conditions for incoming components
	double *Um1 = malloc(m * sizeof(double));		// U^n_{-1}
	double *UN1 = malloc(m * sizeof(double));		// U^n_{N+1}
	for (i = 0; i < m/2; i++)
	{
		// outgoing at right boundary
		UN1[i] = Un[i+m*(N-1)];
		// incoming at left boundary (Maxwell)
		Um1[i] = (1-lambda)*Un[m-1-i] + lambda*fluxL*UmaxwL[i];
	}
	for (i = m/2; i < m; i++)
	{
		// outgoing at left boundary
		Um1[i] = Un[i];
		// incoming at right boundary (Maxwell)
		UN1[i] = (1-lambda)*Un[m-1-i+m*(N-1)] + lambda*fluxR*UmaxwR[i];
	}

	double *Vn     = malloc((N+4) * sizeof(double));
	double *alpha  = malloc((N+3) * sizeof(double));
	double *h_beta = malloc((N+2) * sizeof(double));
	// shift index zero to allow for negative indices
	Vn += 2;
	alpha += 2;
	h_beta++;

	for (i = 0; i < m; i++)
	{
		const double Aabs = fabs(A[i]);

		// nu_p = dt*lambda_p/h
		const double nu = dt*A[i]/h;

		// for slope-limiter term
		const double S    = 0.5 * A[i] * (copysign(1, nu) - nu);
		const double Sabs = 0.5 * Aabs * (copysign(1, nu) - nu);

		Vn[-2] = Um1[i];
		Vn[-1] = Um1[i];
		// local copy
		for (j = 0; j < (int)N; j++) {
			Vn[j] = Un[i+m*j];
		}
		Vn[N]   = UN1[i];
		Vn[N+1] = UN1[i];

		// Eq. (16.34)
		// alpha_j = R^{-1}*(U^n_{j+1} - U^n_j)
		for (j = -2; j <= (int)N; j++) {
			alpha[j] = Vn[j+1] - Vn[j];
		}

		// h * beta_j, Eq. (16.56)
		for (j = -1; j <= (int)N; j++) {
			h_beta[j] = minmod(alpha[j], alpha[j-1]);
		}

		// first term in square brackets of Eq. (16.57) is basically Godunov step in Eq. (13.17)
		const double scale = dt/(2*h);
		for (j = 0; j < (int)N; j++)
		{
			Un1[i+m*j] = Vn[j] + scale * (
				- A[i]*(    Vn[j+1] -     Vn[j-1]) + Aabs*(    Vn[j+1] - 2*    Vn[j] +     Vn[j-1])
				- S   *(h_beta[j+1] - h_beta[j-1]) + Sabs*(h_beta[j+1] - 2*h_beta[j] + h_beta[j-1]));
		}
	}

	// clean up
	h_beta--;
	alpha -= 2;
	Vn -= 2;
	free(h_beta);
	free(alpha);
	free(Vn);
	free(UN1);
	free(Um1);
}
