/// \file simulation.h
/// \brief Header file for simulating the quantum Boltzmann equation, using time splitting to deal with collision, convection and external magnetic field separately
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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "collision.h"

//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data for collision time step
///
typedef struct
{
	CcInterm_t intermCc;	//!< data for Cc collision operator
	CdInterm_t intermCd;	//!< data for Cd collision operator
	wignerF_t *Cc;			//!< conservative collision operator
	wignerF_t *Cd;			//!< dissipative collision operator
	wignerF_t *CB;			//!< external magnetic field
	wignerF_t *Wmid;		//!< intermediate Wigner state for midpoint rule
}
collisionInterm_t;


void CollisionInterm_Create(collisionInterm_t *interm);

void CollisionInterm_Delete(collisionInterm_t *interm);


void CollisionStep(const wignerF_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double dt, const double Bext[3], collisionInterm_t *restrict_ interm, wignerF_t *restrict_ Wnext);


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data for inhomogeneous simulation time step
///
typedef struct
{
	unsigned int numVol;			//!< number of finite volumes
	collisionInterm_t intermColl;	//!< intermediate data for collision step
	double vgrid[N_GRID];			//!< discretized velocities
	grid_t vxgridMaxw;				//!< discretized x-direction velocity grid for transport with Maxwell boundary condition
	wignerV_t *Wtmp;				//!< intermediate Wigner state after first half of transport term
	wignerF_t *WtmpF;				//!< corresponding Fourier representation
	wignerF_t *WtmpF2;				//!< intermediate Wigner state after collision step
	cgrid_t   *WFh;					//!< intermediate state for Fourier transform
	fftw_plan *plan_forw;			//!< array of forward Fourier transform plans
	fftw_plan *plan_back;			//!< array of backward Fourier transform plans
}
inhomStepInterm_t;


void InhomStepInterm_Create(const unsigned int numVol, const double L, inhomStepInterm_t *interm);

void InhomStepInterm_Delete(inhomStepInterm_t *interm);


void InhomStepPeriodic(const wignerV_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wnext);


void InhomStepDirichlet(const wignerV_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wnext);


void InhomStepMaxwell(const wignerV_t *restrict_ W, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double lambda, const wignerV_t *restrict_ WMaxwL, const wignerV_t *restrict_ WMaxwR, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wnext);


//_______________________________________________________________________________________________________________________
///
/// \brief Intermediate data for spatially homogeneous simulation
///
typedef struct
{
	collisionInterm_t intermColl;	//!< intermediate data for collision step
	grid_t    *Wtmp;				//!< temporary state for Fourier transform
	cgrid_t   *WFh;					//!< intermediate state for Fourier transform
	wignerF_t *WtmpF;				//!< actual Fourier representation
	wignerF_t *WtmpF2;				//!< intermediate Wigner state after collision step
	fftw_plan plan_forw;			//!< forward Fourier transform plan
	fftw_plan plan_back[4];			//!< backward Fourier transform plans
}
simulationHomInterm_t;


void SimulationHomInterm_Create(simulationHomInterm_t *interm);

void SimulationHomInterm_Delete(simulationHomInterm_t *interm);


void SimulationHomogeneous(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double dt, const unsigned int numsteps, const double Bext[3], simulationHomInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv);


//_______________________________________________________________________________________________________________________
//


void SimulationInhomogeneousPeriodic(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const unsigned int numsteps, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv);

void SimulationInhomogeneousDirichlet(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const unsigned int numsteps, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv);

void SimulationInhomogeneousMaxwell(const wignerV_t *restrict_ W0, const quadI1_t *restrict_ quadI1, const quadI2_t *restrict_ quadI2, const quadI3_t *restrict_ quadI3, const quadI1_t *restrict_ quadI4,
	const double h, const double dt, const double lambda, const wignerV_t *restrict_ WMaxwL, const wignerV_t *restrict_ WMaxwR, const unsigned int numsteps, const double *restrict_ Bext, inhomStepInterm_t *restrict_ interm, wignerV_t *restrict_ Wevolv);



#endif
