/// \file simulation.c
/// \brief Simulation of the quantum Boltzmann equation using time splitting to deal with collision, convection and external magnetic field separately
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

#include "simulation.h"
#include "finite_volume.h"
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <assert.h>
#ifdef USE_MPI
#include <mpi.h>
#endif


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data for collision time step
///
void CollisionInterm_Create(collisionInterm_t *interm)
{
	CcInterm_Create(&interm->intermCc);
	CdInterm_Create(&interm->intermCd);
	interm->Cc   = fftw_malloc(sizeof(interm->Cc[0]));
	interm->Cd   = fftw_malloc(sizeof(interm->Cd[0]));
	interm->CB   = fftw_malloc(sizeof(interm->CB[0]));
	interm->Wmid = fftw_malloc(sizeof(interm->Wmid[0]));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data for collision time step
///
void CollisionInterm_Delete(collisionInterm_t *interm)
{
	fftw_free(interm->Wmid);
	fftw_free(interm->CB);
	fftw_free(interm->Cd);
	fftw_free(interm->Cc);
	interm->Wmid = NULL;
	interm->Cd   = NULL;
	interm->Cc   = NULL;

	CdInterm_Delete(&interm->intermCd);
	CcInterm_Delete(&interm->intermCc);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Rotation induced by external magnetic field
///
static inline
void MagneticRotation(const wignerF_t *restrict_ W, const double Bext[3], wignerF_t *restrict_ CB)
{
	unsigned int k;

	// set trace terms to zero
	memset(CB, 0, sizeof(CB[0]));

	for (k = 0; k < N_GRID*N_GRID; k++)
	{
		// -i [ B, W ] = 2 ( B cross \vec w_1 )\cdot\sigma

		// y, z
		CB->comp[1].data[k][0] = 2*(Bext[1]*W->comp[3].data[k][0] - Bext[2]*W->comp[2].data[k][0]);
		CB->comp[1].data[k][1] = 2*(Bext[1]*W->comp[3].data[k][1] - Bext[2]*W->comp[2].data[k][1]);

		// z, x
		CB->comp[2].data[k][0] = 2*(Bext[2]*W->comp[1].data[k][0] - Bext[0]*W->comp[3].data[k][0]);
		CB->comp[2].data[k][1] = 2*(Bext[2]*W->comp[1].data[k][1] - Bext[0]*W->comp[3].data[k][1]);

		// x, y
		CB->comp[3].data[k][0] = 2*(Bext[0]*W->comp[2].data[k][0] - Bext[1]*W->comp[1].data[k][0]);
		CB->comp[3].data[k][1] = 2*(Bext[0]*W->comp[2].data[k][1] - Bext[1]*W->comp[1].data[k][1]);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Collision time step using midpoint rule
///
void CollisionStep(const wignerF_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double dt, const double Bext[3], collisionInterm_t *restrict_ interm, wignerF_t *restrict_ Wnext)
{
	unsigned int j, k;

	// midpoint rule

	CcInt(W, quadI1, &interm->intermCc, interm->Cc);
	CdInt(W, quadI2, quadI3, quadI4, &interm->intermCd, interm->Cd);
	MagneticRotation(W, Bext, interm->CB);
	const double dth = 0.5 * dt;
	for (j = 0; j < 4; j++)
	{
		for (k = 0; k < N_GRID*N_GRID; k++)
		{
			interm->Wmid->comp[j].data[k][0] = W->comp[j].data[k][0] + dth * (interm->Cc->comp[j].data[k][0] + interm->Cd->comp[j].data[k][0] + interm->CB->comp[j].data[k][0]);
			interm->Wmid->comp[j].data[k][1] = W->comp[j].data[k][1] + dth * (interm->Cc->comp[j].data[k][1] + interm->Cd->comp[j].data[k][1] + interm->CB->comp[j].data[k][1]);
		}
	}

	CcInt(interm->Wmid, quadI1, &interm->intermCc, interm->Cc);
	CdInt(interm->Wmid, quadI2, quadI3, quadI4, &interm->intermCd, interm->Cd);
	MagneticRotation(interm->Wmid, Bext, interm->CB);
	for (j = 0; j < 4; j++)
	{
		for (k = 0; k < N_GRID*N_GRID; k++)
		{
			Wnext->comp[j].data[k][0] = W->comp[j].data[k][0] + dt * (interm->Cc->comp[j].data[k][0] + interm->Cd->comp[j].data[k][0] + interm->CB->comp[j].data[k][0]);
			Wnext->comp[j].data[k][1] = W->comp[j].data[k][1] + dt * (interm->Cc->comp[j].data[k][1] + interm->Cd->comp[j].data[k][1] + interm->CB->comp[j].data[k][1]);
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data for simulation time step
///
void InhomStepInterm_Create(const unsigned int numVol, const double L, inhomStepInterm_t *interm)
{
	unsigned int i, j, k, l;

	CollisionInterm_Create(&interm->intermColl);

	interm->numVol = numVol;

	const double scale = (2*L)/N_GRID;
	for (k = 0; k < N_GRID; k++) {
		interm->vgrid[k] = k < N_GRID/2 ? scale*k : scale*((int)k - N_GRID);
	}
	// Maxwell velocity grid: omit velocity zero and add velocity +L
	for (l = 0; l < N_GRID; l++) {
		for (k = 0; k < N_GRID; k++) {
			// x components depend on slow index 'l' for compatibility with Matlab's "meshgrid"
			interm->vxgridMaxw.data[k + N_GRID*l] = l < N_GRID/2 ? scale*(l+1) : scale*((int)l - N_GRID);
		}
	}

	interm->Wtmp   = fftw_malloc(numVol * sizeof(interm->Wtmp[0]));
	interm->WtmpF  = fftw_malloc(numVol * sizeof(interm->WtmpF[0]));
	interm->WtmpF2 = fftw_malloc(numVol * sizeof(interm->WtmpF2[0]));
	interm->WFh    = fftw_malloc(sizeof(interm->WFh[0]));

	// create Fourier plans
	interm->plan_forw = malloc(numVol*4 * sizeof(fftw_plan));
	interm->plan_back = malloc(numVol*4 * sizeof(fftw_plan));
	for (i = 0; i < numVol; i++)
	{
		for (j = 0; j < 4; j++)		// for each component...
		{
			interm->plan_forw[j+4*i] = fftw_plan_dft_r2c_2d(N_GRID, N_GRID, interm->Wtmp[i].comp[j].data, interm->WFh->data, FFTW_ESTIMATE);
			interm->plan_back[j+4*i] = fftw_plan_dft_2d(N_GRID, N_GRID, interm->WtmpF2[i].comp[j].data, interm->WtmpF[i].comp[j].data, FFTW_BACKWARD, FFTW_ESTIMATE);
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data for simulation time step
///
void InhomStepInterm_Delete(inhomStepInterm_t *interm)
{
	unsigned int i, j;

	for (i = 0; i < interm->numVol; i++)
	{
		for (j = 0; j < 4; j++)
		{
			fftw_destroy_plan(interm->plan_back[j+4*i]);
			fftw_destroy_plan(interm->plan_forw[j+4*i]);
		}
	}
	free(interm->plan_back);
	free(interm->plan_forw);

	fftw_free(interm->WFh);
	fftw_free(interm->WtmpF2);
	fftw_free(interm->WtmpF);
	fftw_free(interm->Wtmp);

	CollisionInterm_Delete(&interm->intermColl);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Internal collision step function, including forward and backward Fourier transform
///
static inline
void InternalFourierCollisionStep(const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double dt, const double *restrict_ Bext, const bool skipL, const bool skipR, inhomStepInterm_t *restrict_ interm)
{
	unsigned int i, j, k, l;

	const unsigned int i0 = (skipL ? 1 : 0);
	const unsigned int i1 = (skipR ? interm->numVol - 1 : interm->numVol);

	// represent in Fourier space
	for (i = i0; i < i1; i++)
	{
		for (j = 0; j < 4; j++)		// for each component...
		{
			// FFT of real data fills half the resulting complex array only
			// fftw_plan_dft_r2c_2d(N_GRID, N_GRID, interm->Wtmp[i].comp[j].data, interm->WFh->data, 0);
			fftw_execute(interm->plan_forw[j+4*i]);

			const double scale = 1.0 / (N_GRID * N_GRID);
			for (k = 0; k < N_GRID*(N_GRID/2+1); k++)
			{
				interm->WFh->data[k][0] *= scale;
				interm->WFh->data[k][1] *= scale;
			}

			// copy entries to fill complete complex array
			// FFTW uses row-major ordering
			for (l = 0; l < N_GRID; l++)
			{
				// complementary index
				unsigned int lc = (N_GRID - l) & N_MASK;

				for (k = 0; k <= N_GRID/2; k++)
				{
					interm->WtmpF[i].comp[j].data[N_GRID*l + k][0] = interm->WFh->data[(N_GRID/2+1)*l + k][0];
					interm->WtmpF[i].comp[j].data[N_GRID*l + k][1] = interm->WFh->data[(N_GRID/2+1)*l + k][1];
				}
				for (k = N_GRID/2+1; k < N_GRID; k++)
				{
					// complementary index
					unsigned int kc = N_GRID - k;
					// copy conjugate value
					interm->WtmpF[i].comp[j].data[N_GRID*l + k][0] =  interm->WFh->data[(N_GRID/2+1)*lc + kc][0];
					interm->WtmpF[i].comp[j].data[N_GRID*l + k][1] = -interm->WFh->data[(N_GRID/2+1)*lc + kc][1];
				}
			}
		}
	}

	// local collisions (and magnetic rotation) at each x
	for (i = i0; i < i1; i++)
	{
		CollisionStep(&interm->WtmpF[i], quadI1, quadI2, quadI3, quadI4, dt, &Bext[3*i], &interm->intermColl, &interm->WtmpF2[i]);
	}

	// transform to physical velocity space (backward Fourier transform)
	for (i = i0; i < i1; i++)
	{
		for (j = 0; j < 4; j++)		// for each component...
		{
			// fftw_plan_dft_2d(N_GRID, N_GRID, interm->WtmpF2[i].comp[j].data, interm->WtmpF[i].comp[j].data, FFTW_BACKWARD, 0);
			fftw_execute(interm->plan_back[j+4*i]);

			// take physical real part
			for (k = 0; k < N_GRID*N_GRID; k++)
			{
				interm->Wtmp[i].comp[j].data[k] = interm->WtmpF[i].comp[j].data[k][0];
			}
		}
	}
}


//_______________________________________________________________________________________________________________________
//


#ifdef USE_MPI


#define TAG_TRANSPORT_TO_LEFT_H1		1		//!< MPI tag for transferring Wigner states to the left  during first  half of transport calculation
#define TAG_TRANSPORT_TO_LEFT_H2		2		//!< MPI tag for transferring Wigner states to the left  during second half of transport calculation
#define TAG_TRANSPORT_TO_RIGHT_H1		3		//!< MPI tag for transferring Wigner states to the right during first  half of transport calculation
#define TAG_TRANSPORT_TO_RIGHT_H2		4		//!< MPI tag for transferring Wigner states to the right during second half of transport calculation


#endif	// USE_MPI


//_______________________________________________________________________________________________________________________
///
/// \brief Inhomogeneous simulation time step, one-dimensional in space with periodic boundary conditions
///
void InhomStepPeriodic(const wignerV_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wnext)
{
	unsigned int i, j, k, l;

	#ifdef USE_MPI

	// each MPI node handles a subinterval of the spatial domain, and 'W' is assumed to point to exactly this (discretized) interval

	// need to copy two finite volume cells to each neighbor for transport term
	assert(interm->numVol >= 2);

	// +4 for Wigner states from left and right MPI neighbors
	double *Wx  = fftw_malloc((interm->numVol + 4) * sizeof(double));
	double *Wx1 = fftw_malloc((interm->numVol + 4) * sizeof(double));

	int hr;
	MPI_Status status;
	MPI_Request req_left, req_right;

	// size of the computing group
	int groupsize;
	MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
	// current "rank"
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// assume that neighbors have 'rank - 1' and 'rank + 1'
	const int rank_left  = (rank - 1 + groupsize) % groupsize;
	const int rank_right = (rank + 1            ) % groupsize;

	wignerV_t *Wleft  = fftw_malloc(2 * sizeof(wignerV_t));
	wignerV_t *Wright = fftw_malloc(2 * sizeof(wignerV_t));

	// send Wigner states at border to neighboring MPI computing nodes
	hr = MPI_Isend((void *)W, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_LEFT_H1, MPI_COMM_WORLD, &req_left);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Isend()' for sending Wigner state at left border to neighbor with MPI rank %d failed, exiting...\n", rank_left);
		exit(-1);
	}
	hr = MPI_Isend((void *)&W[interm->numVol-2], 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_RIGHT_H1, MPI_COMM_WORLD, &req_right);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Isend()' for sending Wigner state at right border to neighbor with MPI rank %d failed, exiting...\n", rank_right);
		exit(-1);
	}

	// obtain neighboring Wigner states
	// left state arrives via transport to the right
	hr = MPI_Recv(Wleft, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_RIGHT_H1, MPI_COMM_WORLD, &status);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_left);
		exit(-1);
	}
	//// check how many numbers were actually received
	//int count;
	//MPI_Get_count(&status, MPI_DOUBLE, &count);
	//if (count != 2*sizeof(wignerV_t) / sizeof(double))
	//{
	//	fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d has wrong count %d, exiting...\n", rank_left, count);
	//	exit(-1);
	//}
	// right state arrives via transport to the left
	hr = MPI_Recv(Wright, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_LEFT_H1, MPI_COMM_WORLD, &status);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_right);
		exit(-1);
	}
	//// check how many numbers were actually received
	//MPI_Get_count(&status, MPI_DOUBLE, &count);
	//if (count != 2*sizeof(wignerV_t) / sizeof(double))
	//{
	//	fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d has wrong count %d, exiting...\n", rank_right, count);
	//	exit(-1);
	//}

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i + 2] = W[i].comp[j].data[k + N_GRID*l];
				}
				// data from neighbors
				Wx[0]                = Wleft[0] .comp[j].data[k + N_GRID*l];
				Wx[1]                = Wleft[1] .comp[j].data[k + N_GRID*l];
				Wx[interm->numVol+2] = Wright[0].comp[j].data[k + N_GRID*l];
				Wx[interm->numVol+3] = Wright[1].comp[j].data[k + N_GRID*l];

				// transport with velocity v_x for half a time step
				SlopeLimiterStepPeriodic(h, 0.5*dt, interm->vgrid[l], interm->numVol + 4, Wx, Wx1);

				// copy values to temporary variable
				for (i = 0; i < interm->numVol; i++)
				{
					interm->Wtmp[i].comp[j].data[k + N_GRID*l] = Wx1[i + 2];
				}
			}
		}
	}

	// by now the transfer should be completed
	MPI_Wait(&req_left,  &status);
	MPI_Wait(&req_right, &status);

	#else	// !defined(USE_MPI)

	double *Wx  = fftw_malloc(interm->numVol * sizeof(double));
	double *Wx1 = fftw_malloc(interm->numVol * sizeof(double));

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i] = W[i].comp[j].data[k + N_GRID*l];
				}

				// transport with velocity v_x for half a time step
				SlopeLimiterStepPeriodic(h, 0.5*dt, interm->vgrid[l], interm->numVol, Wx, Wx1);

				// copy values to temporary variable
				for (i = 0; i < interm->numVol; i++)
				{
					interm->Wtmp[i].comp[j].data[k + N_GRID*l] = Wx1[i];
				}
			}
		}
	}

	#endif	// USE_MPI

	// local collisions, including forward and backward Fourier transform
	InternalFourierCollisionStep(quadI1, quadI2, quadI3, quadI4, dt, Bext, false, false, interm);

	#ifdef USE_MPI

	hr = MPI_Isend((void *)interm->Wtmp, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_LEFT_H2, MPI_COMM_WORLD, &req_left);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Isend()' for sending Wigner state at left border to neighbor with MPI rank %d failed, exiting...\n", rank_left);
		exit(-1);
	}
	hr = MPI_Isend((void *)&interm->Wtmp[interm->numVol-2], 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_RIGHT_H2, MPI_COMM_WORLD, &req_right);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Isend()' for sending Wigner state at right border to neighbor with MPI rank %d failed, exiting...\n", rank_right);
		exit(-1);
	}
	
	// obtain neighboring Wigner states
	// left state arrives via transport to the right
	hr = MPI_Recv(Wleft, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_RIGHT_H2, MPI_COMM_WORLD, &status);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_left);
		exit(-1);
	}
	// right state arrives via transport to the left
	hr = MPI_Recv(Wright, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_LEFT_H2, MPI_COMM_WORLD, &status);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_right);
		exit(-1);
	}

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i + 2] = interm->Wtmp[i].comp[j].data[k + N_GRID*l];
				}
				// data from neighbors
				Wx[0]                = Wleft[0] .comp[j].data[k + N_GRID*l];
				Wx[1]                = Wleft[1] .comp[j].data[k + N_GRID*l];
				Wx[interm->numVol+2] = Wright[0].comp[j].data[k + N_GRID*l];
				Wx[interm->numVol+3] = Wright[1].comp[j].data[k + N_GRID*l];

				// transport with velocity v_x for half a time step
				SlopeLimiterStepPeriodic(h, 0.5*dt, interm->vgrid[l], interm->numVol + 4, Wx, Wx1);

				// copy values to output Wigner state
				for (i = 0; i < interm->numVol; i++)
				{
					Wnext[i].comp[j].data[k + N_GRID*l] = Wx1[i + 2];
				}
			}
		}
	}

	// by now the transfer should be completed
	MPI_Wait(&req_left,  &status);
	MPI_Wait(&req_right, &status);

	// clean up
	fftw_free(Wright);
	fftw_free(Wleft);

	#else	// !defined(USE_MPI)

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i] = interm->Wtmp[i].comp[j].data[k + N_GRID*l];
				}

				// transport with velocity v_x for half a time step
				SlopeLimiterStepPeriodic(h, 0.5*dt, interm->vgrid[l], interm->numVol, Wx, Wx1);

				// copy values to output Wigner state
				for (i = 0; i < interm->numVol; i++)
				{
					Wnext[i].comp[j].data[k + N_GRID*l] = Wx1[i];
				}
			}
		}
	}

	#endif	// USE_MPI

	fftw_free(Wx1);
	fftw_free(Wx);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Inhomogeneous simulation time step, one-dimensional in space with Dirichlet boundary conditions
///
void InhomStepDirichlet(const wignerV_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wnext)
{
	unsigned int i, j, k, l;

	#ifdef USE_MPI

	// each MPI node handles a subinterval of the spatial domain, and 'W' is assumed to point to exactly this (discretized) interval

	// need to copy two finite volume cells to each neighbor for transport term
	assert(interm->numVol >= 2);

	int hr;
	MPI_Status status;
	MPI_Request req_left, req_right;

	// size of the computing group
	int groupsize;
	MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
	// current "rank"
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// assume that neighbors have 'rank - 1' and 'rank + 1'
	const int rank_left  = rank > 0           ? rank - 1 : rank;
	const int rank_right = rank < groupsize-1 ? rank + 1 : rank;

	// effective number of finite volumes for current computing node;
	// +4 for Wigner states from left and right MPI neighbors
	int numVolEff = interm->numVol + 4;
	if (rank_left  == rank) numVolEff -= 2;		// no left neighbor at left border
	if (rank_right == rank) numVolEff -= 2;		// no right neighbor at right border

	double *Wx  = fftw_malloc(numVolEff * sizeof(double));
	double *Wx1 = fftw_malloc(numVolEff * sizeof(double));

	wignerV_t *Wleft  = fftw_malloc(2 * sizeof(wignerV_t));
	wignerV_t *Wright = fftw_malloc(2 * sizeof(wignerV_t));

	// send Wigner states at border to neighboring MPI computing nodes
	if (rank_left < rank) {
		hr = MPI_Isend((void *)W, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_LEFT_H1, MPI_COMM_WORLD, &req_left);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at left border to neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	if (rank_right > rank) {
		hr = MPI_Isend((void *)&W[interm->numVol-2], 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_RIGHT_H1, MPI_COMM_WORLD, &req_right);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at right border to neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}

	// obtain neighboring Wigner states
	// left state arrives via transport to the right
	if (rank_left < rank) {
		hr = MPI_Recv(Wleft, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_RIGHT_H1, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	// right state arrives via transport to the left
	if (rank_right > rank) {
		hr = MPI_Recv(Wright, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_LEFT_H1, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}

	// half step of transport term (Trotter splitting)
	const int leftShift = rank_left < rank ? 2 : 0;
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i + leftShift] = W[i].comp[j].data[k + N_GRID*l];
				}
				// data from neighbors
				if (rank_left < rank) {
					Wx[0] = Wleft[0].comp[j].data[k + N_GRID*l];
					Wx[1] = Wleft[1].comp[j].data[k + N_GRID*l];
				}
				if (rank_right > rank) {
					Wx[interm->numVol+leftShift  ] = Wright[0].comp[j].data[k + N_GRID*l];
					Wx[interm->numVol+leftShift+1] = Wright[1].comp[j].data[k + N_GRID*l];
				}

				// transport with velocity v_x for half a time step
				SlopeLimiterStepDirichlet(h, 0.5*dt, interm->vgrid[l], numVolEff, Wx, Wx1);

				// copy values to temporary variable
				for (i = 0; i < interm->numVol; i++)
				{
					interm->Wtmp[i].comp[j].data[k + N_GRID*l] = Wx1[i + leftShift];
				}
			}
		}
	}

	// by now the transfer should be completed
	if (rank_left  < rank) MPI_Wait(&req_left,  &status);
	if (rank_right > rank) MPI_Wait(&req_right, &status);

	#else	// !defined(USE_MPI)

	double *Wx  = fftw_malloc(interm->numVol * sizeof(double));
	double *Wx1 = fftw_malloc(interm->numVol * sizeof(double));

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i] = W[i].comp[j].data[k + N_GRID*l];
				}

				// transport with velocity v_x for half a time step
				SlopeLimiterStepDirichlet(h, 0.5*dt, interm->vgrid[l], interm->numVol, Wx, Wx1);

				// copy values to temporary variable
				for (i = 0; i < interm->numVol; i++)
				{
					interm->Wtmp[i].comp[j].data[k + N_GRID*l] = Wx1[i];
				}
			}
		}
	}

	#endif	// USE_MPI

	// local collisions, including forward and backward Fourier transform;
	// omit finite volumes on boundary from collision step
	#ifdef USE_MPI
	InternalFourierCollisionStep(quadI1, quadI2, quadI3, quadI4, dt, Bext, rank_left == rank, rank_right == rank, interm);
	#else
	InternalFourierCollisionStep(quadI1, quadI2, quadI3, quadI4, dt, Bext, true, true, interm);
	#endif

	#ifdef USE_MPI

	if (rank_left < rank) {
		hr = MPI_Isend((void *)interm->Wtmp, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_LEFT_H2, MPI_COMM_WORLD, &req_left);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at left border to neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	if (rank_right > rank) {
		hr = MPI_Isend((void *)&interm->Wtmp[interm->numVol-2], 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_RIGHT_H2, MPI_COMM_WORLD, &req_right);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at right border to neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}
	
	// obtain neighboring Wigner states
	// left state arrives via transport to the right
	if (rank_left < rank) {
		hr = MPI_Recv(Wleft, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_RIGHT_H2, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	// right state arrives via transport to the left
	if (rank_right > rank) {
		hr = MPI_Recv(Wright, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_LEFT_H2, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i + leftShift] = interm->Wtmp[i].comp[j].data[k + N_GRID*l];
				}
				// data from neighbors
				if (rank_left < rank) {
					Wx[0] = Wleft[0].comp[j].data[k + N_GRID*l];
					Wx[1] = Wleft[1].comp[j].data[k + N_GRID*l];
				}
				if (rank_right > rank) {
					Wx[interm->numVol+leftShift  ] = Wright[0].comp[j].data[k + N_GRID*l];
					Wx[interm->numVol+leftShift+1] = Wright[1].comp[j].data[k + N_GRID*l];
				}

				// transport with velocity v_x for half a time step
				SlopeLimiterStepDirichlet(h, 0.5*dt, interm->vgrid[l], numVolEff, Wx, Wx1);

				// copy values to output Wigner state
				for (i = 0; i < interm->numVol; i++)
				{
					Wnext[i].comp[j].data[k + N_GRID*l] = Wx1[i + leftShift];
				}
			}
		}
	}

	// by now the transfer should be completed
	if (rank_left  < rank) MPI_Wait(&req_left,  &status);
	if (rank_right > rank) MPI_Wait(&req_right, &status);

	// clean up
	fftw_free(Wright);
	fftw_free(Wleft);

	#else	// !defined(USE_MPI)

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		for (l = 0; l < N_GRID; l++)
		{
			for (k = 0; k < N_GRID; k++)
			{
				// copy x values into array
				for (i = 0; i < interm->numVol; i++)
				{
					Wx[i] = interm->Wtmp[i].comp[j].data[k + N_GRID*l];
				}

				// transport with velocity v_x for half a time step
				SlopeLimiterStepDirichlet(h, 0.5*dt, interm->vgrid[l], interm->numVol, Wx, Wx1);

				// copy values to output Wigner state
				for (i = 0; i < interm->numVol; i++)
				{
					Wnext[i].comp[j].data[k + N_GRID*l] = Wx1[i];
				}
			}
		}
	}

	#endif	// USE_MPI

	fftw_free(Wx1);
	fftw_free(Wx);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Copy values defined on standard velocity grid to Maxwell velocity grid
///
static inline void CopyToMaxwellVelocity(const grid_t *restrict_ comp, grid_t *restrict_ compMaxw)
{
	// skip components with zero velocity in x-direction;
	// entry with velocity -L in x direction is copied twice to -L and +L (periodicity with period 2*L)
	memcpy( compMaxw->data,                  &comp->data[N_GRID],          sizeof(grid_t)/2);
	memcpy(&compMaxw->data[N_GRID*N_GRID/2], &comp->data[N_GRID*N_GRID/2], sizeof(grid_t)/2);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Copy values defined on Maxwell velocity grid to standard velocity grid
///
static inline void CopyFromMaxwellVelocity(const grid_t *restrict_ compMaxw, const double *restrict_ comp_v0, grid_t *restrict_ comp)
{
	// components with zero velocity in x-direction are stored in 'comp_v0'
	memcpy( comp->data,                   comp_v0,                         N_GRID*sizeof(double));
	memcpy(&comp->data[N_GRID],           compMaxw->data,                  sizeof(grid_t)/2 - N_GRID*sizeof(double));
	memcpy(&comp->data[N_GRID*N_GRID/2], &compMaxw->data[N_GRID*N_GRID/2], sizeof(grid_t)/2);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Transform values defined on standard velocity grid to Maxwell velocity grid
///
static inline void TransformToMaxwellVelocity(grid_t *restrict_ comp)
{
	unsigned int l;
	for (l = 0; l < N_GRID/2; l++)
	{
		// overwrite components with zero velocity in x-direction;
		// entry with velocity -L in x direction is copied twice to -L and +L (periodicity with period 2*L)
		memcpy(&comp->data[l*N_GRID], &comp->data[(l + 1)*N_GRID], N_GRID*sizeof(double));
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Inhomogeneous simulation time step, one-dimensional in space with Maxwell boundary conditions
///
void InhomStepMaxwell(const wignerV_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double lambda, const wignerV_t *restrict_ WMaxwL, const wignerV_t *restrict_ WMaxwR, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wnext)
{
	unsigned int i, j;

	// modified Wigner states for Maxwell reflection operator:
	// omit components with velocity zero in x-direction, and introduce +L velocity components containing a copy of -L
	wignerV_t WMaxwLeff, WMaxwReff;
	for (j = 0; j < 4; j++)
	{
		CopyToMaxwellVelocity(&WMaxwL->comp[j], &WMaxwLeff.comp[j]);
		CopyToMaxwellVelocity(&WMaxwR->comp[j], &WMaxwReff.comp[j]);
	}
	// normalization: only incoming trace component relevant
	{
		double normL = 1.0 / CalculateRightwardFlux(N_GRID*N_GRID, interm->vxgridMaxw.data, WMaxwLeff.comp[0].data);	// rightward flux at left boundary
		double normR = 1.0 / CalculateLeftwardFlux (N_GRID*N_GRID, interm->vxgridMaxw.data, WMaxwReff.comp[0].data);	// leftward flux at right boundary

		// multiply with normalising factor
		for (j = 0; j < 4; j++)
		{
			for (i = 0; i < N_GRID*N_GRID; i++)
			{
				WMaxwLeff.comp[j].data[i] *= normL;
				WMaxwReff.comp[j].data[i] *= normR;
			}
		}
	}

	// calculate outgoing flux based on trace component
	double fluxL, fluxR;
	{
		grid_t WtrL, WtrR;
		CopyToMaxwellVelocity(&W[0               ].comp[0], &WtrL);		// trace component at left  boundary
		CopyToMaxwellVelocity(&W[interm->numVol-1].comp[0], &WtrR);		// trace component at right boundary
		fluxL = CalculateLeftwardFlux (N_GRID*N_GRID, interm->vxgridMaxw.data, WtrL.data);
		fluxR = CalculateRightwardFlux(N_GRID*N_GRID, interm->vxgridMaxw.data, WtrR.data);
	}

	#ifdef USE_MPI

	// each MPI node handles a subinterval of the spatial domain, and 'W' is assumed to point to exactly this (discretized) interval

	// need to copy two finite volume cells to each neighbor for transport term
	assert(interm->numVol >= 2);

	int hr;
	MPI_Status status;
	MPI_Request req_left, req_right;

	// size of the computing group
	int groupsize;
	MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
	// current "rank"
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// assume that neighbors have 'rank - 1' and 'rank + 1'
	const int rank_left  = rank > 0           ? rank - 1 : rank;
	const int rank_right = rank < groupsize-1 ? rank + 1 : rank;

	// effective number of finite volumes for current computing node;
	// +4 for Wigner states from left and right MPI neighbors
	int numVolEff = interm->numVol + 4;
	if (rank_left  == rank) { numVolEff -= 2; }		// no left neighbor at left border
	if (rank_right == rank) { numVolEff -= 2; }		// no right neighbor at right border

	grid_t *Wx  = fftw_malloc(numVolEff * sizeof(grid_t));
	grid_t *Wx1 = fftw_malloc(numVolEff * sizeof(grid_t));

	wignerV_t *Wleft  = fftw_malloc(2 * sizeof(wignerV_t));
	wignerV_t *Wright = fftw_malloc(2 * sizeof(wignerV_t));

	// send Wigner states at border to neighboring MPI computing nodes
	if (rank_left < rank) {
		hr = MPI_Isend((void *)W, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_LEFT_H1, MPI_COMM_WORLD, &req_left);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at left border to neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	if (rank_right > rank) {
		hr = MPI_Isend((void *)&W[interm->numVol-2], 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_RIGHT_H1, MPI_COMM_WORLD, &req_right);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at right border to neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}

	// obtain neighboring Wigner states
	// left state arrives via transport to the right
	if (rank_left < rank) {
		hr = MPI_Recv(Wleft, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_RIGHT_H1, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	// right state arrives via transport to the left
	if (rank_right > rank) {
		hr = MPI_Recv(Wright, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_LEFT_H1, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}

	// half step of transport term (Trotter splitting)
	const int leftShift = rank_left < rank ? 2 : 0;
	for (j = 0; j < 4; j++)
	{
		// copy x values into array
		for (i = 0; i < interm->numVol; i++)
		{
			CopyToMaxwellVelocity(&W[i].comp[j], &Wx[i + leftShift]);
		}
		// data from neighbors
		if (rank_left < rank) {
			CopyToMaxwellVelocity(&Wleft[0].comp[j], &Wx[0]);
			CopyToMaxwellVelocity(&Wleft[1].comp[j], &Wx[1]);
		}
		if (rank_right > rank) {
			CopyToMaxwellVelocity(&Wright[0].comp[j], &Wx[interm->numVol+leftShift  ]);
			CopyToMaxwellVelocity(&Wright[1].comp[j], &Wx[interm->numVol+leftShift+1]);
		}

		// transport for half a time step; intermediate subintervals not affected by Maxwell boundary conditions due to padded states from neighbors
		SlopeLimiterStepMaxwell(h, 0.5*dt, N_GRID*N_GRID, interm->vxgridMaxw.data, numVolEff, lambda, fluxL, fluxR, WMaxwLeff.comp[j].data, WMaxwReff.comp[j].data, Wx->data, Wx1->data);

		// copy values to temporary variable;
		// take components with velocity zero in x-direction from original input Wigner state 'W'
		for (i = 0; i < interm->numVol; i++)
		{
			CopyFromMaxwellVelocity(&Wx1[i + leftShift], W[i].comp[j].data, &interm->Wtmp[i].comp[j]);
		}
	}

	// by now the transfer should be completed
	if (rank_left  < rank) MPI_Wait(&req_left,  &status);
	if (rank_right > rank) MPI_Wait(&req_right, &status);

	#else	// !defined(USE_MPI)

	grid_t *Wx  = fftw_malloc(interm->numVol * sizeof(grid_t));
	grid_t *Wx1 = fftw_malloc(interm->numVol * sizeof(grid_t));

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		// copy x values into array
		for (i = 0; i < interm->numVol; i++)
		{
			CopyToMaxwellVelocity(&W[i].comp[j], &Wx[i]);
		}

		// transport for half a time step
		SlopeLimiterStepMaxwell(h, 0.5*dt, N_GRID*N_GRID, interm->vxgridMaxw.data, interm->numVol, lambda, fluxL, fluxR, WMaxwLeff.comp[j].data, WMaxwReff.comp[j].data, Wx->data, Wx1->data);

		// copy values to temporary variable;
		// take components with velocity zero in x-direction from original input Wigner state 'W'
		for (i = 0; i < interm->numVol; i++)
		{
			CopyFromMaxwellVelocity(&Wx1[i], W[i].comp[j].data, &interm->Wtmp[i].comp[j]);
		}
	}

	#endif	// USE_MPI

	// local collisions, including forward and backward Fourier transform
	InternalFourierCollisionStep(quadI1, quadI2, quadI3, quadI4, dt, Bext, false, false, interm);

	// calculate new outgoing flux based on trace component
	{
		grid_t WtrL, WtrR;
		CopyToMaxwellVelocity(&interm->Wtmp[0               ].comp[0], &WtrL);		// trace component at left  boundary
		CopyToMaxwellVelocity(&interm->Wtmp[interm->numVol-1].comp[0], &WtrR);		// trace component at right boundary
		fluxL = CalculateLeftwardFlux (N_GRID*N_GRID, interm->vxgridMaxw.data, WtrL.data);
		fluxR = CalculateRightwardFlux(N_GRID*N_GRID, interm->vxgridMaxw.data, WtrR.data);
	}

	#ifdef USE_MPI

	if (rank_left < rank) {
		hr = MPI_Isend((void *)interm->Wtmp, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_LEFT_H2, MPI_COMM_WORLD, &req_left);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at left border to neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	if (rank_right > rank) {
		hr = MPI_Isend((void *)&interm->Wtmp[interm->numVol-2], 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_RIGHT_H2, MPI_COMM_WORLD, &req_right);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Isend()' for sending Wigner state at right border to neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}
	
	// obtain neighboring Wigner states
	// left state arrives via transport to the right
	if (rank_left < rank) {
		hr = MPI_Recv(Wleft, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_left, TAG_TRANSPORT_TO_RIGHT_H2, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining left Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_left);
			exit(-1);
		}
	}
	// right state arrives via transport to the left
	if (rank_right > rank) {
		hr = MPI_Recv(Wright, 2*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, rank_right, TAG_TRANSPORT_TO_LEFT_H2, MPI_COMM_WORLD, &status);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Recv()' for obtaining right Wigner state from neighbor with MPI rank %d failed, exiting...\n", rank_right);
			exit(-1);
		}
	}

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		// copy x values into array
		for (i = 0; i < interm->numVol; i++)
		{
			CopyToMaxwellVelocity(&interm->Wtmp[i].comp[j], &Wx[i + leftShift]);
		}
		// data from neighbors
		if (rank_left < rank) {
			CopyToMaxwellVelocity(&Wleft[0].comp[j], &Wx[0]);
			CopyToMaxwellVelocity(&Wleft[1].comp[j], &Wx[1]);
		}
		if (rank_right > rank) {
			CopyToMaxwellVelocity(&Wright[0].comp[j], &Wx[interm->numVol+leftShift  ]);
			CopyToMaxwellVelocity(&Wright[1].comp[j], &Wx[interm->numVol+leftShift+1]);
		}

		// transport for half a time step; intermediate subintervals not affected by Maxwell boundary conditions due to padded states from neighbors
		SlopeLimiterStepMaxwell(h, 0.5*dt, N_GRID*N_GRID, interm->vxgridMaxw.data, numVolEff, lambda, fluxL, fluxR, WMaxwLeff.comp[j].data, WMaxwReff.comp[j].data, Wx->data, Wx1->data);

		// copy values to output Wigner state;
		// take components with velocity zero in x-direction from Wigner state 'interm->Wtmp[i]' before transport
		for (i = 0; i < interm->numVol; i++)
		{
			CopyFromMaxwellVelocity(&Wx1[i + leftShift], interm->Wtmp[i].comp[j].data, &Wnext[i].comp[j]);
		}
	}

	// by now the transfer should be completed
	if (rank_left  < rank) MPI_Wait(&req_left,  &status);
	if (rank_right > rank) MPI_Wait(&req_right, &status);

	// clean up
	fftw_free(Wright);
	fftw_free(Wleft);

	#else	// !defined(USE_MPI)

	// half step of transport term (Trotter splitting)
	for (j = 0; j < 4; j++)
	{
		// copy x values into array
		for (i = 0; i < interm->numVol; i++)
		{
			CopyToMaxwellVelocity(&interm->Wtmp[i].comp[j], &Wx[i]);
		}

		// transport for half a time step
		SlopeLimiterStepMaxwell(h, 0.5*dt, N_GRID*N_GRID, interm->vxgridMaxw.data, interm->numVol, lambda, fluxL, fluxR, WMaxwLeff.comp[j].data, WMaxwReff.comp[j].data, Wx->data, Wx1->data);

		// copy values to output Wigner state;
		// take components with velocity zero in x-direction from Wigner state 'interm->Wtmp[i]' before transport
		for (i = 0; i < interm->numVol; i++)
		{
			CopyFromMaxwellVelocity(&Wx1[i], interm->Wtmp[i].comp[j].data, &Wnext[i].comp[j]);
		}
	}

	#endif	// USE_MPI

	fftw_free(Wx1);
	fftw_free(Wx);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate intermediate data for spatially homogeneous simulation
///
void SimulationHomInterm_Create(simulationHomInterm_t *interm)
{
	unsigned int j;

	CollisionInterm_Create(&interm->intermColl);

	interm->Wtmp   = fftw_malloc(sizeof(interm->Wtmp[0]));
	interm->WFh    = fftw_malloc(sizeof(interm->WFh[0]));
	interm->WtmpF  = fftw_malloc(sizeof(interm->WtmpF[0]));
	interm->WtmpF2 = fftw_malloc(sizeof(interm->WtmpF2[0]));

	// create Fourier plans
	interm->plan_forw = fftw_plan_dft_r2c_2d(N_GRID, N_GRID, interm->Wtmp->data, interm->WFh->data, FFTW_ESTIMATE);
	for (j = 0; j < 4; j++)
	{
		interm->plan_back[j] = fftw_plan_dft_2d(N_GRID, N_GRID, interm->WtmpF2->comp[j].data, interm->WtmpF->comp[j].data, FFTW_BACKWARD, FFTW_ESTIMATE);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of intermediate data for spatially homogeneous simulation
///
void SimulationHomInterm_Delete(simulationHomInterm_t *interm)
{
	unsigned int j;

	for (j = 0; j < 4; j++)
	{
		fftw_destroy_plan(interm->plan_back[j]);
	}
	fftw_destroy_plan(interm->plan_forw);

	fftw_free(interm->WtmpF2);
	fftw_free(interm->WtmpF);
	fftw_free(interm->WFh);
	fftw_free(interm->Wtmp);

	CollisionInterm_Delete(&interm->intermColl);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Spatially homogeneous simulation, store time evolution in output structure 'Wevolv'
///
void SimulationHomogeneous(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double dt, const unsigned int numsteps, const double Bext[3], simulationHomInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv)
{
	unsigned int it, j, k, l;

	// copy initial state
	memcpy(Wevolv, W0, sizeof(wignerV_t));

	for (it = 1; it < numsteps; it++)
	{
		// represent in Fourier space
		for (j = 0; j < 4; j++)		// for each component...
		{
			// FFT of real data fills half the resulting complex array only
			memcpy(interm->Wtmp, &Wevolv[it-1].comp[j].data, sizeof(interm->Wtmp[0]));
			fftw_execute(interm->plan_forw);

			const double scale = 1.0 / (N_GRID * N_GRID);
			for (k = 0; k < N_GRID*(N_GRID/2+1); k++)
			{
				interm->WFh->data[k][0] *= scale;
				interm->WFh->data[k][1] *= scale;
			}

			// copy entries to fill complete complex array
			// FFTW uses row-major ordering
			for (l = 0; l < N_GRID; l++)
			{
				// complementary index
				unsigned int lc = (N_GRID - l) & N_MASK;

				for (k = 0; k <= N_GRID/2; k++)
				{
					interm->WtmpF->comp[j].data[N_GRID*l + k][0] = interm->WFh->data[(N_GRID/2+1)*l + k][0];
					interm->WtmpF->comp[j].data[N_GRID*l + k][1] = interm->WFh->data[(N_GRID/2+1)*l + k][1];
				}
				for (k = N_GRID/2+1; k < N_GRID; k++)
				{
					// complementary index
					unsigned int kc = N_GRID - k;
					// copy conjugate value
					interm->WtmpF->comp[j].data[N_GRID*l + k][0] =  interm->WFh->data[(N_GRID/2+1)*lc + kc][0];
					interm->WtmpF->comp[j].data[N_GRID*l + k][1] = -interm->WFh->data[(N_GRID/2+1)*lc + kc][1];
				}
			}
		}

		// actual collision step in Fourier space
		CollisionStep(interm->WtmpF, quadI1, quadI2, quadI3, quadI4, dt, Bext, &interm->intermColl, interm->WtmpF2);

		// transform to physical velocity space (backward Fourier transform)
		for (j = 0; j < 4; j++)		// for each component...
		{
			fftw_execute(interm->plan_back[j]);

			// take physical real part
			for (k = 0; k < N_GRID*N_GRID; k++)
			{
				Wevolv[it].comp[j].data[k] = interm->WtmpF->comp[j].data[k][0];
			}
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Spatially inhomogeneous simulation with periodic boundary conditions, store time evolution in output structure 'Wevolv'
///
void SimulationInhomogeneousPeriodic(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const unsigned int numsteps, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv)
{
	unsigned int it;

	// copy initial state
	memcpy(Wevolv, W0, interm->numVol * sizeof(wignerV_t));

	for (it = 1; it < numsteps; it++)
	{
		// actual time step
		InhomStepPeriodic(&Wevolv[(it-1)*interm->numVol], quadI1, quadI2, quadI3, quadI4, h, dt, Bext, interm, &Wevolv[it*interm->numVol]);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Spatially inhomogeneous simulation with Dirichlet boundary conditions, store time evolution in output structure 'Wevolv'
///
void SimulationInhomogeneousDirichlet(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const unsigned int numsteps, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv)
{
	unsigned int it;

	// copy initial state
	memcpy(Wevolv, W0, interm->numVol * sizeof(wignerV_t));

	for (it = 1; it < numsteps; it++)
	{
		// actual time step
		InhomStepDirichlet(&Wevolv[(it-1)*interm->numVol], quadI1, quadI2, quadI3, quadI4, h, dt, Bext, interm, &Wevolv[it*interm->numVol]);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Spatially inhomogeneous simulation with Maxwell boundary conditions, store time evolution in output structure 'Wevolv'
///
void SimulationInhomogeneousMaxwell(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double lambda, const wignerV_t *restrict_ WMaxwL, const wignerV_t *restrict_ WMaxwR, const unsigned int numsteps, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv)
{
	unsigned int it;

	// copy initial state
	memcpy(Wevolv, W0, interm->numVol * sizeof(wignerV_t));

	for (it = 1; it < numsteps; it++)
	{
		// actual time step
		InhomStepMaxwell(&Wevolv[(it-1)*interm->numVol], quadI1, quadI2, quadI3, quadI4, h, dt, lambda, WMaxwL, WMaxwR, Bext, interm, &Wevolv[it*interm->numVol]);
	}
}
