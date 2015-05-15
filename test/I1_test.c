// Test file for calculating the I1 integral
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
	const unsigned int J = 32;
	const double L = 12;
	const double R = 7.5;

	unsigned int i;
	int hr;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	quadI1_t quad;
	FourierI1(J, L, R, &quad);

	cgrid_t *f = fftw_malloc(sizeof(cgrid_t));
	cgrid_t *g = fftw_malloc(sizeof(cgrid_t));
	cgrid_t *h = fftw_malloc(sizeof(cgrid_t));

	I1interm_t interm;
	I1interm_Create(&interm);

	cgrid_t *I1     = fftw_malloc(sizeof(cgrid_t));
	cgrid_t *I1_ref = fftw_malloc(sizeof(cgrid_t));

	// load data from disk

	V_RETURN(ReadData("../test/data/f.dat", f->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));
	V_RETURN(ReadData("../test/data/g.dat", g->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));
	V_RETURN(ReadData("../test/data/h.dat", h->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));

	printf("calculating I1 integral...\n");

	// perform computation
	I1integral(f, g, h, &quad, &interm, I1);

	printf("done.\n\n");

	// example
	printf("example:\n");
	printf("I1->data[0]: %g + %gi\n", I1->data[0][0], I1->data[0][1]);
	printf("I1->data[1]: %g + %gi\n", I1->data[1][0], I1->data[1][1]);

	// compare with reference
	V_RETURN(ReadData("../test/data/I1_ref.dat", I1_ref->data, sizeof(fftw_complex), sizeof(cgrid_t) / sizeof(fftw_complex)));
	double err = 0;
	for (i = 0; i < sizeof(I1_ref->data) / sizeof(fftw_complex); i++)
	{
		err += fabs(I1->data[i][0] - I1_ref->data[i][0]) + fabs(I1->data[i][1] - I1_ref->data[i][1]);
	}
	printf("\naverage error: %g\n", err / (sizeof(I1_ref->data) / sizeof(fftw_complex)));

	// clean up

	I1interm_Delete(&interm);

	fftw_free(I1_ref);
	fftw_free(I1);
	fftw_free(h);
	fftw_free(g);
	fftw_free(f);

	QuadI1_Delete(&quad);

	fftw_cleanup();

	return 0;
}
