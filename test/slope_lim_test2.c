// Test file for simulating a finite volume time step using slope limiters, with Maxwell boundary conditions
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
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }



int main()
{
	const double h  = 0.01;		// spatial mesh width
	const double dt = 0.005;	// time step

	const unsigned int m = 4;
	const unsigned int N = 101;

	const double A[4] = { 1.0, 0.5, -0.5, -1.0 };

	// accommodation coefficient for Maxwell reflection operator
	const double lambda = 0.4;

	// predetermined fluxes
	const double fluxL = 0.7;
	const double fluxR = 0.45;

	unsigned int j;
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// incoming states at the left and right boundary for diffusive reflection operator
	double UmaxwL[4] = { 0.2, -M_PI_4, 0,        0   };		// non-zero for velocities > 0
	double UmaxwR[4] = { 0,    0,      1.0/M_E, -0.8 };		// non-zero for velocities < 0
	//// normalization
	//double normL = 1.0 / CalculateRightwardFlux(m, A, UmaxwL);
	//double normR = 1.0 / CalculateLeftwardFlux (m, A, UmaxwR);
	//for (j = 0; j < m; j++)
	//{
	//	UmaxwL[j] *= normL;
	//	UmaxwR[j] *= normR;
	//}

	// load data from disk
	double *Un = malloc(m*N * sizeof(double));
	V_RETURN(ReadData("../test/data/UnM.dat", Un, sizeof(double), m*N));

	double *Un1     = malloc(m*N * sizeof(double));
	double *Un1_ref = malloc(m*N * sizeof(double));

	printf("simulating a finite volume time step using slope limiters, with Maxwell boundary conditions...\n");

	// slope limiter step
	SlopeLimiterStepMaxwell(h, dt, m, A, N, lambda, fluxL, fluxR, UmaxwL, UmaxwR, Un, Un1);

	printf("done.\n\n");

	// example
	printf("example:\n");
	printf("Un1[0]:  %g\n", Un1[0]);
	printf("Un1[17]: %g\n", Un1[17]);

	// load reference data from disk
	V_RETURN(ReadData("../test/data/Un1_Maxwell_ref.dat", Un1_ref, sizeof(double), m*N));
	printf("reference Un1[0]:  %g\n", Un1_ref[0]);
	printf("reference Un1[17]: %g\n", Un1_ref[17]);

	// compare with reference
	double err = 0;
	for (j = 0; j < m*N; j++)
	{
		err += fabs(Un1[j] - Un1_ref[j]);
	}
	printf("\ncumulative error: %g\n", err);

	// clean up
	free(Un1_ref);
	free(Un1);
	free(Un);

	return 0;
}
