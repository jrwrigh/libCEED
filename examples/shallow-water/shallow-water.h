// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

/// @file
/// Initial condition and operator for the shallow-water equations example using PETSc

#ifndef densitycurrent_h
#define densitycurrent_h

#ifndef CeedPragmaOMP
#  ifdef _OPENMP
#    define CeedPragmaOMP_(a) _Pragma(#a)
#    define CeedPragmaOMP(a) CeedPragmaOMP_(omp a)
#  else
#    define CeedPragmaOMP(a)
#  endif
#endif

#include <math.h>

// *****************************************************************************
// This QFunction sets the the initial conditions and boundary conditions
//
//  TO DO
//
// *****************************************************************************
static int ICsSW(void *ctx, CeedInt Q,
                 const CeedScalar *const *in, CeedScalar *const *out) {

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

  // Inputs
  const CeedScalar *X = in[0];
  // Outputs
  CeedScalar *q0 = out[0], *coords = out[1];
  // Context
  const CeedScalar *context = (const CeedScalar*)ctx;
  const CeedScalar u0     = context[0];
  const CeedScalar v0     = context[1];
  const CeedScalar h0     = context[2];

  CeedPragmaOMP(simd)
  // Quadrature Point Loop
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    // -- Coordinates
    const CeedScalar x = X[i+0*Q];
    const CeedScalar y = X[i+1*Q];

    // Initial Conditions
    q0[i+0*Q] = u0;
    q0[i+1*Q] = v0;
    q0[i+2*Q] = h0;

    // Coordinates
    coords[i+0*Q] = x;
    coords[i+1*Q] = y;

  } // End of Quadrature Point Loop

  // Return
  return 0;
}

// *****************************************************************************
// This QFunction implements the explicit terms of the shallow-water
// equations
//
// The equations represent 2D shallow-water flow on a spherical surface, where
// the state variable, h, represents the height function.
//
// State (scalar) variable: h
//
// Shallow-water Equations explicit terms:
//   TO DO
// *****************************************************************************
static int SWExplicit(void *ctx, CeedInt Q,
              const CeedScalar *const *in, CeedScalar *const *out) {
  // Inputs
  const CeedScalar *q = in[0], *dq = in[1], *qdata = in[2], *x = in[3];
  // Outputs
  CeedScalar *v = out[0], *dv = out[1];
  // Context
  const CeedScalar *context = (const CeedScalar*)ctx;

  CeedPragmaOMP(simd)
  // Quadrature Point Loop
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    // -- Interp in
    const CeedScalar u[2]   = { q[i+0*Q] / rho,
                                q[i+1*Q] / rho
                                };
    const CeedScalar h      =   q[i+2*Q];
        // -- Grad in
    const CeedScalar du[4]  =  { dq[i+(0+3*0)*Q],
                                 dq[i+(0+3*1)*Q],
                                 dq[i+(1+3*0)*Q],
                                 dq[i+(2+3*1)*Q]
                                };
    const CeedScalar dh[2]  = { dq[i+(2+3*0)*Q],
                                dq[i+(2+3*1)*Q]
                               };
    // -- Interp-to-Interp qdata
    const CeedScalar wJ     =   qdata[i+ 0*Q];
    // -- Interp-to-Grad qdata
    //      Symmetric 3x3 matrix
    const CeedScalar wBJ[4] = { qdata[i+ 1*Q],
                                qdata[i+ 2*Q],
                                qdata[i+ 3*Q],
                                qdata[i+ 4*Q]
                              };
    // -- Grad-to-Grad qdata
    const CeedScalar wBBJ[3]  = { qdata[i+5*Q],
                                  qdata[i+6*Q],
                                  qdata[i+7*Q]
                                };

    // -- curl u
    const CeedScalar curlu[2] = { du[2],
                                 -du[1]
                                 };

    // The Physics

    // -- Equation for h
    // ---- u h
    dv[i+0*Q]  = rho*u[0]*wBJ[0] + rho*u[1]*wBJ[1] + rho*u[2]*wBJ[2];
    dv[i+1*Q]  = rho*u[0]*wBJ[3] + rho*u[1]*wBJ[4] + rho*u[2]*wBJ[5];
    dv[i+2*Q]  = rho*u[0]*wBJ[6] + rho*u[1]*wBJ[7] + rho*u[2]*wBJ[8];
    // ---- No Change
    v[i+0*Q] = 0;


  } // End Quadrature Point Loop

  // Return
  return 0;
}

// *****************************************************************************
// This QFunction implements the implicit terms of the shallow-water
// equations
//
// The equations represent 2D shallow-water flow on a spherical surface, where
// the state variable, h, represents the height function.
//
// State (scalar) variable: h
//
// Shallow-water Equations implicit terms:
//   TO DO
// *****************************************************************************
static int SWImplicit(void *ctx, CeedInt Q,
              const CeedScalar *const *in, CeedScalar *const *out) {
  // Inputs
  const CeedScalar *q = in[0], *dq = in[1], *qdata = in[2], *x = in[3];
  // Outputs
  CeedScalar *v = out[0], *dv = out[1];
  // Context
  const CeedScalar *context = (const CeedScalar*)ctx;

  CeedPragmaOMP(simd)
  // Quadrature Point Loop
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    // -- Interp in
    const CeedScalar u[2]   = { q[i+0*Q] / rho,
                                q[i+1*Q] / rho
                                };
    const CeedScalar h      =   q[i+2*Q];
        // -- Grad in
    const CeedScalar du[4]  =  { dq[i+(0+3*0)*Q],
                                 dq[i+(0+3*1)*Q],
                                 dq[i+(1+3*0)*Q],
                                 dq[i+(2+3*1)*Q]
                                };
    const CeedScalar dh[2]  = { dq[i+(2+3*0)*Q],
                                dq[i+(2+3*1)*Q]
                               };
    // -- Interp-to-Interp qdata
    const CeedScalar wJ     =   qdata[i+ 0*Q];
    // -- Interp-to-Grad qdata
    //      Symmetric 3x3 matrix
    const CeedScalar wBJ[4] = { qdata[i+ 1*Q],
                                qdata[i+ 2*Q],
                                qdata[i+ 3*Q],
                                qdata[i+ 4*Q]
                              };
    // -- Grad-to-Grad qdata
    const CeedScalar wBBJ[3]  = { qdata[i+5*Q],
                                  qdata[i+6*Q],
                                  qdata[i+7*Q]
                                };

    // -- curl u
    const CeedScalar curlu[2] = { du[2],
                                 -du[1]
                                 };

    // The Physics

    // -- Equation for h
    // ---- u h
    dv[i+0*Q]  = rho*u[0]*wBJ[0] + rho*u[1]*wBJ[1] + rho*u[2]*wBJ[2];
    dv[i+1*Q]  = rho*u[0]*wBJ[3] + rho*u[1]*wBJ[4] + rho*u[2]*wBJ[5];
    dv[i+2*Q]  = rho*u[0]*wBJ[6] + rho*u[1]*wBJ[7] + rho*u[2]*wBJ[8];
    // ---- No Change
    v[i+0*Q] = 0;


  } // End Quadrature Point Loop

  // Return
  return 0;
}

// *****************************************************************************
#endif
