// Test file for calculating the I2 integral
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

#include "integrals.h"
#include "util.h"


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

#define V_RETURN(x)  { hr = (x); if (hr < 0) { fprintf(stderr, "%s, line %d: command '%s' failed, return value: %d\n", __FILE__, __LINE__, #x, hr); return hr; } }


int main()
{
	const unsigned int J = 7;

	unsigned int i, j;
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	quadI2_t quad;
	QuadI2_Allocate(J, &quad);

	cgrid_t *f = fftw_malloc(sizeof(cgrid_t));
	cgrid_t *g = fftw_malloc(sizeof(cgrid_t));
	cgrid_t *h = fftw_malloc(sizeof(cgrid_t));

	I2interm_t interm;
	I2interm_Create(&interm);

	cgrid_t *I2     = fftw_malloc(sizeof(cgrid_t));
	cgrid_t *I2_ref = fftw_malloc(sizeof(cgrid_t));

	// load data from disk

	V_RETURN(ReadData("../test/data/f.dat", f->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));
	V_RETURN(ReadData("../test/data/g.dat", g->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));
	V_RETURN(ReadData("../test/data/h.dat", h->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));

	//printf("f->data[0]: %g + %gi\n", f->data[0][0], f->data[0][1]);
	//printf("f->data[1]: %g + %gi\n", f->data[1][0], f->data[1][1]);

	//printf("g->data[0]: %g + %gi\n", g->data[0][0], g->data[0][1]);
	//printf("g->data[1]: %g + %gi\n", g->data[1][0], g->data[1][1]);

	//printf("h->data[0]: %g + %gi\n", h->data[0][0], h->data[0][1]);
	//printf("h->data[1]: %g + %gi\n", h->data[1][0], h->data[1][1]);

	double *weight = fftw_malloc(J * sizeof(double));
	V_RETURN(ReadData("../test/data/weight.dat", weight, sizeof(double), J));

	for (j = 0; j < J; j++)
	{
		//printf("weight[%d]: %g\n", j, weight[j]);

		char fname[1024];
		sprintf(fname, "../test/data/psiR1_%d.dat", j + 1);
		V_RETURN(ReadData(fname, &quad.psiR1[j].data, sizeof(double), sizeof(quad.psiR1[j]) / sizeof(double)));
		sprintf(fname, "../test/data/psiR2_%d.dat", j + 1);
		V_RETURN(ReadData(fname, &quad.psiR2[j].data, sizeof(double), sizeof(quad.psiR2[j]) / sizeof(double)));

		// absorb 'weight' into psiR2
		for (i = 0; i < sizeof(quad.psiR2[j]) / sizeof(double); i++)
		{
			quad.psiR2[j].data[i] *= weight[j];
		}
	}

	printf("calculating I2 integral...\n");

	// perform computation
	I2integral(f, g, h, &quad, &interm, I2);

	printf("done.\n\n");

	// example
	printf("example:\n");
	printf("I2->data[0]: %g + %gi\n", I2->data[0][0], I2->data[0][1]);
	printf("I2->data[1]: %g + %gi\n", I2->data[1][0], I2->data[1][1]);

	//WriteData("../test/data/f_psi2F_test.dat", quad.f_psi2F->data, sizeof(fftw_complex), sizeof(quad.f_psi2F->data) / sizeof(fftw_complex), false);

	// compare with reference
	V_RETURN(ReadData("../test/data/I2_ref.dat", I2_ref->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));
	double err = 0;
	for (i = 0; i < sizeof(I2_ref->data) / sizeof(fftw_complex); i++)
	{
		err += fabs(I2->data[i][0] - I2_ref->data[i][0]) + fabs(I2->data[i][1] - I2_ref->data[i][1]);
	}
	printf("\naverage error: %g\n", err / (sizeof(I2_ref->data) / sizeof(fftw_complex)));

	// clean up

	I2interm_Delete(&interm);

	fftw_free(weight);

	fftw_free(I2_ref);
	fftw_free(I2);
	fftw_free(h);
	fftw_free(g);
	fftw_free(f);

	QuadI2_Delete(&quad);

	fftw_cleanup();

	return 0;
}
