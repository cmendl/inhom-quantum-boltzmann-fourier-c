// Test file for simulating a finite volume time step using slope limiters, with periodic and Dirichlet boundary conditions
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
#include <memory.h>
#include <stdio.h>


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }



int main()
{
	const double h  = 0.08;		// spatial mesh width
	const double dt = 0.02;		// time step

	const unsigned int N = 128;

	const double A = -1.2;

	unsigned int j;
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// load data from disk
	double *Un = malloc(N * sizeof(double));
	V_RETURN(ReadData("../test/data/Un.dat", Un, sizeof(double), N));

	double *Un1P     = malloc(N * sizeof(double));
	double *Un1P_ref = malloc(N * sizeof(double));
	double *Un1D     = malloc(N * sizeof(double));
	double *Un1D_ref = malloc(N * sizeof(double));

	printf("simulating a finite volume time step using slope limiters, with periodic and Dirichlet boundary conditions...\n");

	// slope limiter step
	SlopeLimiterStepPeriodic (h, dt, A, N, Un, Un1P);
	SlopeLimiterStepDirichlet(h, dt, A, N, Un, Un1D);

	printf("done.\n\n");

	// example
	printf("example:\n");
	printf("periodic  Un1P[0]:   %g\n", Un1P[0]);
	printf("periodic  Un1P[5]:   %g\n", Un1P[5]);
	printf("periodic  Un1P[N-1]: %g\n", Un1P[N-1]);
	printf("Dirichlet Un1D[0]:   %g\n", Un1D[0]);
	printf("Dirichlet Un1D[5]:   %g\n", Un1D[5]);
	printf("Dirichlet Un1D[N-1]: %g\n", Un1D[N-1]);

	// load reference data from disk
	V_RETURN(ReadData("../test/data/Un1_periodic_ref.dat",  Un1P_ref, sizeof(double), N));
	V_RETURN(ReadData("../test/data/Un1_Dirichlet_ref.dat", Un1D_ref, sizeof(double), N));
	printf("periodic  Un1P_ref[0]:   %g\n", Un1P_ref[0]);
	printf("periodic  Un1P_ref[5]:   %g\n", Un1P_ref[5]);
	printf("periodic  Un1P_ref[N-1]: %g\n", Un1P_ref[N-1]);
	printf("Dirichlet Un1D_ref[0]:   %g\n", Un1D_ref[0]);
	printf("Dirichlet Un1D_ref[5]:   %g\n", Un1D_ref[5]);
	printf("Dirichlet Un1D_ref[N-1]: %g\n", Un1D_ref[N-1]);

	// compare with reference
	double err = 0;
	for (j = 0; j < N; j++)
	{
		err += fabs(Un1P[j] - Un1P_ref[j]) + fabs(Un1D[j] - Un1D_ref[j]);
	}
	printf("\ncumulative error: %g\n", err);

	// clean up
	free(Un1D_ref);
	free(Un1D);
	free(Un1P_ref);
	free(Un1P);
	free(Un);

	return 0;
}
