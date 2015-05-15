/// \file util.h
/// \brief Header file for the utility functions
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

#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>


//_______________________________________________________________________________________________________________________
///
/// \brief square function x -> x^2
///
static inline double square(const double x)
{
	return x*x;
}


//_______________________________________________________________________________________________________________________
///
/// \brief maximum of two numbers
///
static inline double maxf(const double x, const double y)
{
	if (x >= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

//_______________________________________________________________________________________________________________________
///
/// \brief minimum of two numbers
///
static inline double minf(const double x, const double y)
{
	if (x <= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}


//_______________________________________________________________________________________________________________________
//


int ReadData(const char *filename, void *data, const size_t size, const size_t n);

int WriteData(const char *filename, const void *data, const size_t size, const size_t n, const bool append);



#endif
