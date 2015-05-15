/// \file finite_volume.h
/// \brief Header file for the finite volume simulation using slope limiters with various boundary conditions
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

#ifndef FINITE_VOLUME_H
#define FINITE_VOLUME_H


void SlopeLimiterStepPeriodic(const double h, const double dt, const double A, const unsigned int N, const double *restrict_ Un, double *restrict_ Un1);


void SlopeLimiterStepDirichlet(const double h, const double dt, const double A, const unsigned int N, const double *restrict_ Un, double *restrict_ Un1);


double CalculateRightwardFlux(const unsigned int m, const double *restrict_ A, const double *restrict_ U);

double CalculateLeftwardFlux (const unsigned int m, const double *restrict_ A, const double *restrict_ U);

void SlopeLimiterStepMaxwell(const double h, const double dt, const unsigned int m, const double *restrict_ A, const unsigned int N,
	const double lambda, const double fluxL, const double fluxR, const double *restrict_ UmaxwL, const double *restrict_ UmaxwR, const double *restrict_ Un, double *restrict_ Un1);



#endif
