// Test file for simulating a time step of the spatially inhomogeneous Boltzmann equation
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
#include "util.h"
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <assert.h>
#ifdef USE_MPI
#include <mpi.h>
#define TAG_COLLECT			0		//!< MPI tag for collecting results from computing nodes
#endif


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }



#ifdef USE_MPI
int main(int argc, char *argv[])
#else
int main()
#endif
{
	const unsigned int J = 16;
	const unsigned int M = 32;
	const double L = 12;
	const double R = 7.5;
	const unsigned int numVol = 8;

	const double dt = 0.002;	// time step
	const double h  = 0.08;		// spatial mesh width

	unsigned int i, j, k;
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// vanishing external magnetic field
	double *Bext = calloc(3*numVol, sizeof(double));

	// load data from disk
	wignerV_t *W = fftw_malloc(numVol * sizeof(wignerV_t));
	V_RETURN(ReadData("../test/data/Wx.dat", W[0].comp[0].data, sizeof(double), numVol * sizeof(wignerV_t) / sizeof(double)));

	// quadrature rules
	quadI1_t quadI1;
	quadI2_t quadI2;
	quadI3_t quadI3;
	quadI1_t quadI4;
	FourierI1(J, L, R, &quadI1);
	FourierI2(J, L, R, &quadI2);
	FourierI3(J, L, R, M, &quadI3);
	FourierI4(J, L, R, &quadI4);

	// intermediate data for simulation step
	inhomStepInterm_t interm;

	#ifdef USE_MPI

	MPI_Init(&argc, &argv);

	MPI_Status status;

	// size of the computing group
	int groupsize;
	MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
	// current "rank"
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if ((int)numVol < 2*groupsize) {
		fprintf(stderr, "Number of finite volumes (%d) cannot be smaller than twice the number of computing processes (%d), exiting...\n", numVol, groupsize);
		MPI_Finalize();
		exit(-1);
	}

	// distribute spatial finite volumes among computing nodes
	int i0 = numVol *  rank      / groupsize;
	int i1 = numVol * (rank + 1) / groupsize;
	// at least two finite volumnes
	assert(i1 - i0 >= 2);

	// allocate intermediate data
	InhomStepInterm_Create(i1 - i0, L, &interm);

	// next Wigner state on simulation subinterval
	wignerV_t *Wsub_next = fftw_malloc((i1 - i0) * sizeof(wignerV_t));

	printf("simulation time step with %d spatial volumes, process %d...\n", i1 - i0, rank);

	// start timer
	clock_t t_start = clock();

	// inhomogeneous simulation time step
	InhomStepPeriodic(W + i0, &quadI1, &quadI2, &quadI3, &quadI4, h, dt, Bext + 3*i0, &interm, Wsub_next);

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("process %d finished calculation, CPU time: %g\n", rank, cpu_time);

	if (rank == 0)
	{
		wignerV_t *Wnext = fftw_malloc(numVol * sizeof(wignerV_t));

		// collect results from all processes
		memcpy(Wnext, Wsub_next, (i1 - i0) * sizeof(wignerV_t));
		int r;
		for (r = 1; r < groupsize; r++)
		{
			int j0 = numVol *  r      / groupsize;
			int j1 = numVol * (r + 1) / groupsize;
			hr = MPI_Recv(&Wnext[j0], (j1 - j0) * sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, r, TAG_COLLECT, MPI_COMM_WORLD, &status);
			if (hr != MPI_SUCCESS)
			{
				fprintf(stderr, "'MPI_Recv()' for obtaining Wigner states on subinterval from rank %d process failed, exiting...\n", r);
				exit(-1);
			}
		}

		// load reference data from disk
		wignerV_t *Wnext_ref = fftw_malloc(numVol * sizeof(wignerV_t));
		V_RETURN(ReadData("../test/data/Wx_next_ref.dat", Wnext_ref[0].comp[0].data, sizeof(double), numVol * sizeof(wignerV_t) / sizeof(double)));

		// compare with reference
		double err = 0;
		double nrf = 0;
		for (i = 0; i < numVol; i++)
		{
			for (j = 0; j < 4; j++)
			{
				for (k = 0; k < N_GRID*N_GRID; k++)
				{
					err += fabs(Wnext[i].comp[j].data[k] - Wnext_ref[i].comp[j].data[k]);
					nrf += fabs(Wnext_ref[i].comp[j].data[k]);
				}
			}
		}
		printf("cumulative error: %g, relative error: %g\n", err, err / nrf);

		fftw_free(Wnext_ref);
		fftw_free(Wnext);
	}
	else
	{
		// send results to rank 0 process

		hr = MPI_Send(Wsub_next, (i1 - i0)*sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, 0, TAG_COLLECT, MPI_COMM_WORLD);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Send()' for sending Wigner states on subinterval from rank %d to rank 0 process failed, exiting...\n", rank);
			exit(-1);
		}
	}

	fftw_free(Wsub_next);

	#else	// !defined(USE_MPI)

	// allocate intermediate data
	InhomStepInterm_Create(numVol, L, &interm);

	wignerV_t *Wnext = fftw_malloc(numVol * sizeof(wignerV_t));

	printf("simulating a time step of the spatially inhomogeneous Boltzmann equation...\n");

	// start timer
	clock_t t_start = clock();

	// inhomogeneous simulation time step
	InhomStepPeriodic(W, &quadI1, &quadI2, &quadI3, &quadI4, h, dt, Bext, &interm, Wnext);

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("finished calculation, CPU time: %g\n\n", cpu_time);

	// example
	printf("example:\n");
	printf("Wnext[2].comp[0].data[5]: %g\n", Wnext[2].comp[0].data[5]);
	printf("Wnext[2].comp[1].data[5]: %g\n", Wnext[2].comp[1].data[5]);

	// load reference data from disk
	wignerV_t *Wnext_ref = fftw_malloc(numVol * sizeof(wignerV_t));
	V_RETURN(ReadData("../test/data/Wx_next_ref.dat", Wnext_ref[0].comp[0].data, sizeof(double), numVol * sizeof(wignerV_t) / sizeof(double)));
	printf("Wnext_ref[2].comp[0].data[5]: %g\n", Wnext_ref[2].comp[0].data[5]);
	printf("Wnext_ref[2].comp[1].data[5]: %g\n", Wnext_ref[2].comp[1].data[5]);

	// compare with reference
	double err = 0;
	double nrf = 0;
	for (i = 0; i < numVol; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < N_GRID*N_GRID; k++)
			{
				err += fabs(Wnext[i].comp[j].data[k] - Wnext_ref[i].comp[j].data[k]);
				nrf += fabs(Wnext_ref[i].comp[j].data[k]);
			}
		}
	}
	printf("\ncumulative error: %g, relative error: %g\n", err, err / nrf);

	fftw_free(Wnext_ref);
	fftw_free(Wnext);

	#endif

	// clean up
	InhomStepInterm_Delete(&interm);
	QuadI1_Delete(&quadI4);
	QuadI3_Delete(&quadI3);
	QuadI2_Delete(&quadI2);
	QuadI1_Delete(&quadI1);
	fftw_free(W);
	free(Bext);

	fftw_cleanup();

	#ifdef USE_MPI
	MPI_Finalize();
	#endif

	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtDumpMemoryLeaks();
	#endif

	return 0;
}
