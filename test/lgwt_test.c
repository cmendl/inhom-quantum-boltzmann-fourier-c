// Test file for constructing a Legendre-Gauss quadrature rule
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
#include <stdlib.h>
#include <memory.h>


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }


void LegendreGaussQuad(const unsigned int M, const double a, const double b, double *x, double *w);


int main()
{
	const unsigned int M = 15;
	const double a = 0.5;
	const double b = 2;

	unsigned int i;
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	double *x = malloc(M * sizeof(double));
	double *w = malloc(M * sizeof(double));

	printf("constructing a Legendre-Gauss quadrature rule...\n");
	LegendreGaussQuad(M, a, b, x, w);

	// load reference data from disk
	double *x_ref = malloc(M * sizeof(double));
	double *w_ref = malloc(M * sizeof(double));
	V_RETURN(ReadData("../test/data/x_ref.dat", x_ref, sizeof(double), M));
	V_RETURN(ReadData("../test/data/w_ref.dat", w_ref, sizeof(double), M));

	// compare with reference
	double err = 0;
	for (i = 0; i < M; i++)
	{
		err += fabs(x[i] - x_ref[i]) + fabs(w[i] - w_ref[i]);
	}
	printf("total error compared to reference: %g\n", err);

	// clean up

	free(w_ref);
	free(x_ref);

	free(w);
	free(x);

	return 0;
}