/// \file main_inhom.c
/// \brief Main file for simulating the spatially inhomogeneous quantum Boltzmann equation; optional command line parameter is filename of simulation parameter file
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
#include <string.h>
#include <assert.h>
#ifdef USE_MPI
#include <mpi.h>
#define TAG_COLLECT			0		//!< MPI tag for collecting results from computing nodes
#endif


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }


//_______________________________________________________________________________________________________________________
///
/// \brief Boundary types
///
enum boundaryType
{
	BOUNDARY_PERIODIC  = 0,
	BOUNDARY_DIRICHLET = 1,
	BOUNDARY_MAXWELL   = 2,
};


//_______________________________________________________________________________________________________________________
///
/// \brief Simulation parameters
///
typedef struct
{
	unsigned int J;				//!< number of angular discretization points
	unsigned int M;				//!< number of Legendre-Gauss quadrature nodes
	double L;					//!< periodized interval is [-L,L] in each direction
	double R;					//!< Wigner functions should be supported on ball with radius R

	double h;					//!< spatial mesh width
	unsigned int numVol;		//!< number of finite volumes
	double dt;					//!< time step
	unsigned int numsteps;		//!< number of time steps

	char filenameBext[1024];	//!< file name containing external magnetic field; file must store 3*numVol double values (for each spatial finite volume)

	enum boundaryType bType;	//!< boundary type
	double lambda;				//!< accommodation coefficient for reflection operator (only used for Maxwell boundary condition)

	char filenameWinit[1024];	//!< file name of initial Wigner state
	char filenameWMaxwL[1024];	//!< file name of left Wigner state used for diffusive reflection (only for Maxwell boundary condition)
	char filenameWMaxwR[1024];	//!< file name of right Wigner state used for diffusive reflection (only for Maxwell boundary condition)
	char filenameWevolv[1024];	//!< output file name
}
simParams_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Set default parameters
///
static void DefaultParameters(simParams_t *params)
{
	params->J = 32;
	params->M = 32;
	params->L = 12;
	params->R = 7.5;

	params->h = 0.1;	// spatial mesh width
	const double xmax = 1;
	params->numVol = (int)(xmax/params->h);

	params->dt = 0.008;
	const double tmax = 1.0;
	params->numsteps = (int)(tmax / params->dt) + 1;

	// empty file name means zero external magnetic field
	memset(params->filenameBext, 0, sizeof(params->filenameBext));

	params->bType = BOUNDARY_PERIODIC;
	params->lambda = 0.4;	// only used for Maxwell boundary condition

	// default input file name
	strcpy(params->filenameWinit, "../data/W0_inhom.dat");

	memset(params->filenameWMaxwL, 0, sizeof(params->filenameWMaxwL));
	memset(params->filenameWMaxwR, 0, sizeof(params->filenameWMaxwR));

	// default output file name
	strcpy(params->filenameWevolv, "../data/Wevolv_inhom.dat");
}


//_______________________________________________________________________________________________________________________
///
/// \brief Parse parameter text file containing simulation parameters; parameters not set in the file remain unchanged
///
static int ParseParameterFile(const char *filename, simParams_t *params)
{
	FILE *fd = fopen(filename, "r");
	if (fd == NULL) {
		fprintf(stderr, "Cannot open file '%s'.\n", filename);
		return -1;
	}

	// read file line by line
	char line[1024];
	while (fgets(line, sizeof(line), fd) != NULL)
	{
		// check if there is a '='
		if (strchr(line, '=') == NULL) {
			continue;
		}

		// parameter name
		char name[1024];
		strncpy(name, strtok(line, " =\r\n"), 1024);

		// parameter value
		char *value = strtok(NULL, " =\r\n");
		if (value == NULL)
		{
			fprintf(stderr, "Missing value for parameter '%s' in file '%s'.\n", name, filename);
			return -1;
		}

		if (strcmp(name, "J") == 0) {
			params->J = atoi(value);
		}
		else if (strcmp(name, "M") == 0) {
			params->M = atoi(value);
		}
		else if (strcmp(name, "L") == 0) {
			params->L = atoi(value);
		}
		else if (strcmp(name, "R") == 0) {
			params->R = atof(value);
		}
		else if (strcmp(name, "h") == 0) {
			params->h = atof(value);
		}
		else if (strcmp(name, "numVol") == 0) {
			params->numVol = atoi(value);
		}
		else if (strcmp(name, "dt") == 0) {
			params->dt = atof(value);
		}
		else if (strcmp(name, "numsteps") == 0) {
			params->numsteps = atoi(value);
		}
		else if (strcmp(name, "boundaryType") == 0)
		{
			if (strcmp(value, "periodic") == 0) {
				params->bType = BOUNDARY_PERIODIC;
			}
			else if (strcmp(value, "Dirichlet") == 0) {
				params->bType = BOUNDARY_DIRICHLET;
			}
			else if (strcmp(value, "Maxwell") == 0) {
				params->bType = BOUNDARY_MAXWELL;
			}
			else {
				fprintf(stderr, "Warning: unrecognized boundary type '%s' in file '%s'.\n", value, filename);
			}
		}
		else if (strcmp(name, "lambda") == 0) {
			params->lambda= atof(value);
		}
		else if (strcmp(name, "filenameBext") == 0)
		{
			strcpy(params->filenameBext, value);
		}
		else if (strcmp(name, "filenameWinit") == 0)
		{
			strcpy(params->filenameWinit, value);
		}
		else if (strcmp(name, "filenameWMaxwL") == 0)
		{
			strcpy(params->filenameWMaxwL, value);
		}
		else if (strcmp(name, "filenameWMaxwR") == 0)
		{
			strcpy(params->filenameWMaxwR, value);
		}
		else if (strcmp(name, "filenameWevolv") == 0)
		{
			strcpy(params->filenameWevolv, value);
		}
		else {
			fprintf(stderr, "Warning: unrecognized parameter '%s' in file '%s'.\n", name, filename);
		}
	}

	return 0;
}


//_______________________________________________________________________________________________________________________
//


int main(int argc, char *argv[])
{
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// simulation parameters
	simParams_t params;

	#ifdef USE_MPI

	MPI_Init(&argc, &argv);

	// size of the computing group
	int groupsize;
	MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
	// current "rank"
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only rank 0 process reads simulation parameters from file
	if (rank == 0)
	{

	#endif

	// set default parameters
	DefaultParameters(&params);
	// read parameters from file, overwriting default parameters
	if (argc >= 2) {
		printf("Reading simulation parameters from file '%s'...\n", argv[1]);
		hr = ParseParameterFile(argv[1], &params);
		if (hr < 0) {
			fprintf(stderr, "Error parsing parameter file, exiting...\n");
			return -1;
		}
	}

	// print parameters
	printf("Simulation parameters:\n");
	printf("J:           %d\n", params.J);
	printf("M:           %d\n", params.M);
	printf("L:           %g\n", params.L);
	printf("R:           %g\n", params.R);
	printf("h:           %g\n", params.h);
	printf("numVol:      %d\n", params.numVol);
	printf("dt:          %g\n", params.dt);
	printf("numsteps:    %d\n", params.numsteps);
	printf("Bext file:   %s\n", params.filenameBext);
	const char *bTypeString[] = { "periodic", "Dirichlet", "Maxwell" };
	printf("boundary:    %s\n", bTypeString[params.bType]);
	if (params.bType == BOUNDARY_MAXWELL) {
		printf("lambda:      %g\n", params.lambda);
		printf("WMaxwL file: %s\n", params.filenameWMaxwL);
		printf("WMaxwR file: %s\n", params.filenameWMaxwR);
	}
	printf("input file:  %s\n", params.filenameWinit);
	printf("output file: %s\n", params.filenameWevolv);

	#ifdef USE_MPI
	}	// if (rank == 0)

	hr = MPI_Bcast((void *)&params, sizeof(simParams_t), MPI_BYTE, 0, MPI_COMM_WORLD);
	if (hr != MPI_SUCCESS)
	{
		fprintf(stderr, "'MPI_Bcast()' for broadcasting simulation parameters failed, exiting...\n");
		MPI_Finalize();
		return -1;
	}

	#endif

	// external magnetic field
	double *Bext = calloc(3*params.numVol, sizeof(double));
	if (strlen(params.filenameBext) > 0)
	{
		printf("Loading external magnetic field from disk...\n");
		V_RETURN(ReadData(params.filenameBext, Bext, sizeof(double), 3*params.numVol));
	}

	// load data from disk
	printf("Loading initial data from disk...\n");
	wignerV_t *W0 = fftw_malloc(params.numVol * sizeof(wignerV_t));
	V_RETURN(ReadData(params.filenameWinit, W0, sizeof(double), params.numVol * sizeof(wignerV_t) / sizeof(double)));

	// Wigner states for left and right diffusive reflection operator (only used for Maxwell boundary condition)
	wignerV_t WMaxwL, WMaxwR;
	if (params.bType == BOUNDARY_MAXWELL) {
		printf("Loading data for diffusive reflection operator from disk...\n");
		V_RETURN(ReadData(params.filenameWMaxwL, &WMaxwL, sizeof(double), sizeof(wignerV_t) / sizeof(double)));
		V_RETURN(ReadData(params.filenameWMaxwR, &WMaxwR, sizeof(double), sizeof(wignerV_t) / sizeof(double)));
	}

	// quadrature rules
	quadI1_t quadI1;
	quadI2_t quadI2;
	quadI3_t quadI3;
	quadI1_t quadI4;
	FourierI1(params.J, params.L, params.R, &quadI1);
	FourierI2(params.J, params.L, params.R, &quadI2);
	FourierI3(params.J, params.L, params.R, params.M, &quadI3);
	FourierI4(params.J, params.L, params.R, &quadI4);

	// intermediate data
	inhomStepInterm_t interm;

	#ifdef USE_MPI

	MPI_Status status;

	if ((int)params.numVol < 2*groupsize) {
		fprintf(stderr, "Number of finite volumes (%d) cannot be smaller than twice the number of computing processes (%d), exiting...\n", params.numVol, groupsize);
		MPI_Finalize();
		return -1;
	}

	// distribute spatial finite volumes among computing nodes
	int i0 = params.numVol *  rank      / groupsize;
	int nx = params.numVol * (rank + 1) / groupsize - i0;
	// at least two finite volumnes
	assert(nx >= 2);

	// allocate intermediate data
	InhomStepInterm_Create(nx, params.L, &interm);

	// time evolution of Wigner states on simulation subinterval
	wignerV_t *Wsub_evolv = fftw_malloc(params.numsteps*nx * sizeof(wignerV_t));

	printf("Starting simulation with %d spatial volumes, process %d...\n", nx, rank);

	// start timer
	clock_t t_start = clock();

	// spatially inhomogeneous simulation
	if (params.bType == BOUNDARY_PERIODIC) {
		SimulationInhomogeneousPeriodic(W0 + i0, &quadI1, &quadI2, &quadI3, &quadI4, params.h, params.dt, params.numsteps, Bext + 3*i0, &interm, Wsub_evolv);
	}
	else if (params.bType == BOUNDARY_DIRICHLET) {
		SimulationInhomogeneousDirichlet(W0 + i0, &quadI1, &quadI2, &quadI3, &quadI4, params.h, params.dt, params.numsteps, Bext + 3*i0, &interm, Wsub_evolv);
	}
	else if (params.bType == BOUNDARY_MAXWELL) {
		SimulationInhomogeneousMaxwell(W0 + i0, &quadI1, &quadI2, &quadI3, &quadI4, params.h, params.dt, params.lambda, &WMaxwL, &WMaxwR, params.numsteps, Bext + 3*i0, &interm, Wsub_evolv);
	}
	else {
		fprintf(stderr, "Unknown boundary condition, exiting...\n");
		MPI_Finalize();
		return -1;
	}

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("Process %d finished simulation, CPU time: %g\n", rank, cpu_time);

	if (rank == 0)
	{
		wignerV_t *Wevolv = fftw_malloc(params.numsteps*params.numVol * sizeof(wignerV_t));

		unsigned int it;

		// interleave data
		for (it = 0; it < params.numsteps; it++)
		{
			memcpy(&Wevolv[it * params.numVol], &Wsub_evolv[it * nx], nx * sizeof(wignerV_t));
		}

		// collect results from all other processes
		int r;
		for (r = 1; r < groupsize; r++)
		{
			int j0 = params.numVol *  r      / groupsize;
			int ny = params.numVol * (r + 1) / groupsize - j0;

			wignerV_t *Wtmp = fftw_malloc(params.numsteps*ny * sizeof(wignerV_t));

			hr = MPI_Recv(Wtmp, params.numsteps*ny * sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, r, TAG_COLLECT, MPI_COMM_WORLD, &status);
			if (hr != MPI_SUCCESS)
			{
				fprintf(stderr, "'MPI_Recv()' for obtaining Wigner states on subinterval from rank %d process failed, exiting...\n", r);
				MPI_Finalize();
				return -1;
			}

			// interleave data
			for (it = 0; it < params.numsteps; it++)
			{
				memcpy(&Wevolv[it * params.numVol + j0], &Wtmp[it * ny], ny * sizeof(wignerV_t));
			}

			fftw_free(Wtmp);
		}

		// save results to disk
		WriteData(params.filenameWevolv, Wevolv, sizeof(wignerV_t), params.numsteps * params.numVol, false);

		fftw_free(Wevolv);
	}
	else
	{
		// send results to rank 0 process

		hr = MPI_Send(Wsub_evolv, params.numsteps*nx * sizeof(wignerV_t) / sizeof(double), MPI_DOUBLE, 0, TAG_COLLECT, MPI_COMM_WORLD);
		if (hr != MPI_SUCCESS)
		{
			fprintf(stderr, "'MPI_Send()' for sending Wigner states on subinterval from rank %d to rank 0 process failed, exiting...\n", rank);
			MPI_Finalize();
			return -1;
		}
	}

	fftw_free(Wsub_evolv);

	#else	// !defined(USE_MPI)

	// allocate intermediate data
	InhomStepInterm_Create(params.numVol, params.L, &interm);

	// time evolution of Wigner states
	wignerV_t *Wevolv = fftw_malloc(params.numsteps*params.numVol * sizeof(wignerV_t));

	printf("Starting spatially inhomogeneous simulation with %d spatial volumes and %d time steps...\n", params.numVol, params.numsteps);

	// start timer
	clock_t t_start = clock();

	// spatially inhomogeneous simulation
	if (params.bType == BOUNDARY_PERIODIC) {
		SimulationInhomogeneousPeriodic(W0, &quadI1, &quadI2, &quadI3, &quadI4, params.h, params.dt, params.numsteps, Bext, &interm, Wevolv);
	}
	else if (params.bType == BOUNDARY_DIRICHLET) {
		SimulationInhomogeneousDirichlet(W0, &quadI1, &quadI2, &quadI3, &quadI4, params.h, params.dt, params.numsteps, Bext, &interm, Wevolv);
	}
	else if (params.bType == BOUNDARY_MAXWELL) {
		SimulationInhomogeneousMaxwell(W0, &quadI1, &quadI2, &quadI3, &quadI4, params.h, params.dt, params.lambda, &WMaxwL, &WMaxwR, params.numsteps, Bext, &interm, Wevolv);
	}
	else {
		fprintf(stderr, "Unknown boundary condition, exiting...\n");
		return -1;
	}

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("Finished simulation, CPU time: %g\n", cpu_time);

	// save results to disk
	WriteData(params.filenameWevolv, Wevolv, sizeof(wignerV_t), params.numsteps * params.numVol, false);

	#endif

	// clean up
	InhomStepInterm_Delete(&interm);
	QuadI1_Delete(&quadI4);
	QuadI3_Delete(&quadI3);
	QuadI2_Delete(&quadI2);
	QuadI1_Delete(&quadI1);
	fftw_free(W0);
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
