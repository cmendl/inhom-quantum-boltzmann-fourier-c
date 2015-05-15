/// \file main_hom.c
/// \brief Main file for simulating the spatially homogeneous quantum Boltzmann equation
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
#include <time.h>


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }



int main()
{
	const unsigned int J = 32;
	const unsigned int M = 32;
	const double L = 12;
	const double R = 7.5;

	const double dt = 0.001;
	const double tmax = 0.1;
	const unsigned int numsteps = (int)(tmax / dt) + 1;

	// vanishing external magnetic field
	const double Bext[3] = { 0, 0, 0 };

	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// load data from disk
	printf("Loading initial data from disk...\n");
	wignerV_t *W0 = fftw_malloc(sizeof(wignerV_t));
	V_RETURN(ReadData("../data/W0_hom.dat", W0, sizeof(double), sizeof(wignerV_t) / sizeof(double)));

	// quadrature rules
	quadI1_t quadI1;
	quadI2_t quadI2;
	quadI3_t quadI3;
	quadI1_t quadI4;
	FourierI1(J, L, R, &quadI1);
	FourierI2(J, L, R, &quadI2);
	FourierI3(J, L, R, M, &quadI3);
	FourierI4(J, L, R, &quadI4);

	// intermediate data
	simulationHomInterm_t interm;
	SimulationHomInterm_Create(&interm);

	// time evolution of Wigner states
	wignerV_t *Wevolv = fftw_malloc(numsteps * sizeof(wignerV_t));

	printf("Starting spatially homogeneous simulation with %d time steps...\n", numsteps);

	// start timer
	clock_t t_start = clock();

	// homogeneous simulation
	SimulationHomogeneous(W0, &quadI1, &quadI2, &quadI3, &quadI4, dt, numsteps, Bext, &interm, Wevolv);

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("Finished simulation, CPU time: %g\n", cpu_time);

	// save results to disk
	WriteData("../data/Wevolv_hom.dat", Wevolv, sizeof(wignerV_t), numsteps, false);

	// clean up
	fftw_free(Wevolv);
	SimulationHomInterm_Delete(&interm);
	QuadI1_Delete(&quadI4);
	QuadI3_Delete(&quadI3);
	QuadI2_Delete(&quadI2);
	QuadI1_Delete(&quadI1);
	fftw_free(W0);

	fftw_cleanup();

	return 0;
}
