/// \file integrals.h
/// \brief Header file for evaluating the I1 (and I4), I2 and I3 integrals
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

#ifndef INTEGRALS_H
#define INTEGRALS_H

#include "grid.h"
#include "quadrature.h"


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data and FFTW plans for calculating the I1 integral
///
typedef struct
{
	cgrid2_t *f2;			//!< zero-padded f
	cgrid2_t *f2F;			//!< Fourier transform of f2

	cgrid2_t *h2;			//!< zero-padded h
	cgrid2_t *h2F;			//!< Fourier transform of h2

	cgrid_t  *tmp;			//!< temporary variable
	cgrid2_t *tmp2;			//!< temporary variable on extended grid
	cgrid2_t *tmp2F;		//!< Fourier transform of 'tmp2'

	fftw_plan plan_f;		//!< Fourier transform of f2
	fftw_plan plan_h;		//!< Fourier transform of h2
	fftw_plan plan_back;	//!< backward Fourier transform for convolution
	fftw_plan plan_forw;	//!< forward  Fourier transform for convolution
}
I1interm_t;


void I1interm_Create(I1interm_t* interm);

void I1interm_Delete(I1interm_t* interm);


void I1integral(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, const cgrid_t *restrict_ h, const quadI1_t *restrict_ quad, I1interm_t *restrict_ interm, cgrid_t *restrict_ I1);


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data and FFTW plans for calculating the I2 integral
///
typedef struct
{
	cgrid_t *f_psi;			//!< f .* psiR1
	cgrid_t *g_psi;			//!< g .* psiR2

	cgrid2_t *f_psi2;		//!< zero-padded f_psi
	cgrid2_t *g_psi2;		//!< zero-padded g_psi

	cgrid2_t *f_psi2F;		//!< Fourier transform of f_psi2
	cgrid2_t *g_psi2F;		//!< Fourier transform of g_psi2

	cgrid2_t *h2;			//!< zero-padded h
	cgrid2_t *h2F;			//!< Fourier transform of h2

	cgrid_t  *tmp;			//!< temporary variable
	cgrid2_t *tmp2;			//!< temporary variable on extended grid
	cgrid2_t *tmp2F;		//!< Fourier transform of 'tmp2'

	fftw_plan plan_f;		//!< Fourier transform of f_psi2
	fftw_plan plan_g;		//!< Fourier transform of g_psi2
	fftw_plan plan_h;		//!< Fourier transform of h2
	fftw_plan plan_back;	//!< backward Fourier transform for convolution
}
I2interm_t;


void I2interm_Create(I2interm_t* interm);

void I2interm_Delete(I2interm_t* interm);


void I2integral(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, const cgrid_t *restrict_ h, const quadI2_t *restrict_ quad, I2interm_t *restrict_ interm, cgrid_t *restrict_ I2);


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data and FFTW plans for calculating the I3 integral
///
typedef struct
{
	cgrid_t *f_psi;			//!< f .* psiR1
	cgrid_t *h_psi;			//!< h .* psiR2

	cgrid2_t *f_psi2;		//!< zero-padded f_psi
	cgrid2_t *h_psi2;		//!< zero-padded h_psi

	cgrid2_t *f_psi2F;		//!< Fourier transform of f_psi2
	cgrid2_t *h_psi2F;		//!< Fourier transform of h_psi2

	cgrid2_t *g2;			//!< zero-padded g
	cgrid2_t *g2F;			//!< Fourier transform of g2

	cgrid_t  *tmp;			//!< temporary variable
	cgrid2_t *tmp2;			//!< temporary variable on extended grid
	cgrid2_t *tmp2F;		//!< Fourier transform of 'tmp2'

	fftw_plan plan_g;		//!< Fourier transform of g2
	fftw_plan plan_f;		//!< Fourier transform of f_psi2
	fftw_plan plan_h;		//!< Fourier transform of h_psi2
	fftw_plan plan_back;	//!< backward Fourier transform for convolution
	fftw_plan plan_forw;	//!< forward  Fourier transform for convolution
}
I3interm_t;


void I3interm_Create(I3interm_t* interm);

void I3interm_Delete(I3interm_t* interm);


void I3integral(const cgrid_t *restrict_ f, const cgrid_t *restrict_ g, const cgrid_t *restrict_ h, const quadI3_t *restrict_ quad, I3interm_t *restrict_ interm, cgrid_t *restrict_ I3);



#endif
