/// \file collision.h
/// \brief Header file for the collision operator
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

#ifndef COLLISION_H
#define COLLISION_H

#include "integrals.h"
#include "wigner.h"


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data for calculating Cd collision operator
///
typedef struct
{
	I2interm_t intermI2;	//!< intermediate data for I2 integral
	I3interm_t intermI3;	//!< intermediate data for I3 integral
	I1interm_t intermI4;	//!< intermediate data for I4 integral
	wignerF_t *Wt;			//!< 1 - W
	cgrid_t *tmp;			//!< temporary variable
}
CdInterm_t;


void CdInterm_Create(CdInterm_t *interm);

void CdInterm_Delete(CdInterm_t *interm);


void CdInt(const wignerF_t *restrict_ W, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4, CdInterm_t *restrict_ interm, wignerF_t *restrict_ Cd);


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data for calculating Cc collision operator
///
typedef struct
{
	I1interm_t intermI1;	//!< intermediate data for I1 integral
	cgrid_t *Wtr2m1;		//!< 2 W_tr - 1
	cgrid_t *I1;			//!< store I1 integrals
}
CcInterm_t;


void CcInterm_Create(CcInterm_t *interm);

void CcInterm_Delete(CcInterm_t *interm);


void CcInt(const wignerF_t *restrict_ W, const quadI1_t *restrict_ quadI1, CcInterm_t *restrict_ interm, wignerF_t *restrict_ Cc);



#endif
