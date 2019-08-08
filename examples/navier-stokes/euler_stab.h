
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
/// Density current initial condition and operator for Navier-Stokes example using PETSc

// Model from:
//   Semi-Implicit Formulations of the Navier-Stokes Equations: Application to
//   Nonhydrostatic Atmospheric Modeling, Giraldo, Restelli, and Lauter (2010).

#ifndef euler_stab_h
#define euler_stab_h

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
// These initial conditions are given in terms of potential temperature and
//   Exner pressure and then converted to density and total energy.
//   Initial momentum density is zero.
//
// Initial Conditions:
//   Potential Temperature:
//     theta = thetabar + deltatheta
//       thetabar   = theta0 exp( N**2 z / g )
//       deltatheta = r <= rc : theta0(1 + cos(pi r/rc)) / 2
//                     r > rc : 0
//         r        = sqrt( (x - xc)**2 + (y - yc)**2 + (z - zc)**2 )
//         with (xc,yc,zc) center of domain, rc characteristic radius of thermal bubble
//   Exner Pressure:
//     Pi = Pibar + deltaPi
//       Pibar      = g**2 (exp( - N**2 z / g ) - 1) / (cp theta0 N**2)
//       deltaPi    = 0 (hydrostatic balance)
//   Velocity/Momentum Density:
//     Ui = ui = 0
//
// Conversion to Conserved Variables:
//   rho = P0 Pi**(cv/Rd) / (Rd theta)
//   E   = rho (cv theta Pi + (u u)/2 + g z)
//
//  Boundary Conditions:
//    Mass Density:
//      0.0 flux
//    Momentum Density:
//      0.0
//    Energy Density:
//      0.0 flux
//
// Constants:
//   theta0          ,  Potential temperature constant
//   thetaC          ,  Potential temperature perturbation
//   P0              ,  Pressure at the surface
//   N               ,  Brunt-Vaisala frequency
//   cv              ,  Specific heat, constant volume
//   cp              ,  Specific heat, constant pressure
//   Rd     = cp - cv,  Specific heat difference
//   g               ,  Gravity
//   rc              ,  Characteristic radius of thermal bubble
//   lx              ,  Characteristic length scale of domain in x
//   ly              ,  Characteristic length scale of domain in y
//   lz              ,  Characteristic length scale of domain in z
//
// *****************************************************************************
static int ICsEulerStab(void *ctx, CeedInt Q,
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
  const CeedScalar theta0     = context[0];
  const CeedScalar thetaC     = context[1];
  const CeedScalar P0         = context[2];
  const CeedScalar N          = context[3];
  const CeedScalar cv         = context[4];
  const CeedScalar cp         = context[5];
  const CeedScalar Rd         = context[6];
  const CeedScalar g          = context[7];
  const CeedScalar rc         = context[8];
  const CeedScalar lx         = context[9];
  const CeedScalar ly         = context[10];
  const CeedScalar lz         = context[11];
  // Setup
  const CeedScalar tol = 1.e-14;
  const CeedScalar center[3] = {0.5*lx, 0.5*ly, 0.5*lz};

  CeedPragmaOMP(simd)
  // Quadrature Point Loop
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    // -- Coordinates
    const CeedScalar x = X[i+0*Q];
    const CeedScalar y = X[i+1*Q];
    const CeedScalar z = X[i+2*Q];
    // -- Potential temperature, density current
    const CeedScalar r = sqrt(pow((x - center[0]), 2) +
                              pow((y - center[1]), 2) +
                              pow((z - center[2]), 2));
    const CeedScalar deltatheta = r <= rc ? thetaC*(1. + cos(M_PI*r/rc))/2. : 0.;
    const CeedScalar theta = theta0*exp(N*N*z/g) + deltatheta;
    // -- Exner pressure, hydrostatic balance
    const CeedScalar Pi = 1. + g*g*(exp(-N*N*z/g) - 1.) / (cp*theta0*N*N);
    // -- Density
    const CeedScalar rho = P0 * pow(Pi, cv/Rd) / (Rd*theta);

    // Initial Conditions
    q0[i+0*Q] = rho;
    q0[i+1*Q] = 0.0;
    q0[i+2*Q] = 0.0;
    q0[i+3*Q] = 0.0;
    q0[i+4*Q] = rho * (cv*theta*Pi + g*z);

    // Homogeneous Dirichlet Boundary Conditions for Momentum
    if ( fabs(x - 0.0) < tol || fabs(x - lx) < tol ||
         fabs(y - 0.0) < tol || fabs(y - ly) < tol ||
         fabs(z - 0.0) < tol || fabs(z - lz) < tol ) {
      q0[i+1*Q] = 0.0;
      q0[i+2*Q] = 0.0;
      q0[i+3*Q] = 0.0;
    }

    // Coordinates
    coords[i+0*Q] = x;
    coords[i+1*Q] = y;
    coords[i+2*Q] = z;

  } // End of Quadrature Point Loop

  // Return
  return 0;
}

// *****************************************************************************
// This QFunction implements the following formulation of Navier-Stokes
//
// This is 3D compressible Navier-Stokes in conservation form with state
//   variables of density, momentum density, and total energy density.
//
// State Variables: q = ( rho, U1, U2, U3, E )
//   rho - Mass Density
//   Ui  - Momentum Density   ,  Ui = rho ui
//   E   - Total Energy Density,  E  = rho cv T + rho (u u) / 2 + rho g z
//
// Navier-Stokes Equations:
//   drho/dt + div( U )                                = 0
//   dU/dt   + div( rho (u x u) + P I3 ) + rho g khat  = 0
//   dE/dt   + div( (E + P) u )         + rho g u[z]   = 0
//
//
// Equation of State:
//   P = (gamma - 1) (E - rho (u u) / 2)
//********************************************
// Stabilization:
//
//   Tau = [TauC, TauM, TauM, TauM, TauE]
//      f1 = rho  sqrt(2 / (C1  dt) + ui uj gij)    C1 = 1.
//           gij = dXi/dX * dXi/dX 
// TauC = Cc f1 / (8 gii)                           Cc = 1.  
// TauM = 1 / f1
// TauE = TauM / (Ce cv)                            Ce =1.
//
//  SUPG = Galerkin + grad(v) . ( Ai^T * Tau * (Aj q,j) )
//
// Constants:
//   lambda = - 2 / 3,  From Stokes hypothesis
//   mu              ,  Dynamic viscosity
//   k               ,  Thermal conductivity
//   cv              ,  Specific heat, constant volume
//   cp              ,  Specific heat, constant pressure
//   g               ,  Gravity
//   gamma  = cp / cv,  Specific heat ratio: (1 + cp/cv)
//
// *****************************************************************************
static int EulerStab(void *ctx, CeedInt Q,
              const CeedScalar *const *in, CeedScalar *const *out) {
  // Inputs
  const CeedScalar *q = in[0], *dq = in[1], *qdata = in[2], *x = in[3];
  // Outputs
  CeedScalar *v = out[0], *dv = out[1];
  // Context
  const CeedScalar *context = (const CeedScalar*)ctx;
  const CeedScalar k          = context[2];
  const CeedScalar cv         = context[3];
  const CeedScalar cp         = context[4];
  const CeedScalar g          = context[5];
  const CeedScalar Rd         = context[6];
  const CeedScalar gamma      = cp / cv;

  CeedPragmaOMP(simd)
  // Quadrature Point Loop
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    // -- Interp in
    const CeedScalar rho      =    q[i+0*Q];
    const CeedScalar u[3]     =  { q[i+1*Q] / rho,
                                   q[i+2*Q] / rho,
                                   q[i+3*Q] / rho
                                 };
    const CeedScalar E        =    q[i+4*Q];

    // -- Grad in
    const CeedScalar drho[3]  =  { dq[i+(0+5*0)*Q],
                                   dq[i+(0+5*1)*Q],
                                   dq[i+(0+5*2)*Q]
                                 };
    const CeedScalar dU[3][3] = { {dq[i+(1+5*0)*Q],
                                   dq[i+(1+5*1)*Q],
                                   dq[i+(1+5*2)*Q]},
                                  {dq[i+(2+5*0)*Q],
                                   dq[i+(2+5*1)*Q],
                                   dq[i+(2+5*2)*Q]},
                                  {dq[i+(3+5*0)*Q],
                                   dq[i+(3+5*1)*Q],
                                   dq[i+(3+5*2)*Q]}
                                 };
    const CeedScalar dE[3]    =  { dq[i+(4+5*0)*Q],
                                   dq[i+(4+5*1)*Q],
                                   dq[i+(4+5*2)*Q]
                                 };

    // -- Interp-to-Interp qdata
    const CeedScalar wJ       =    qdata[i+ 0*Q];
    // -- Interp-to-Grad qdata
    //      Symmetric 3x3 matrix
    const CeedScalar wBJ[9]   =  { qdata[i+1*Q],
                                   qdata[i+2*Q],
                                   qdata[i+3*Q],
                                   qdata[i+4*Q],
                                   qdata[i+5*Q],
                                   qdata[i+6*Q],
                                   qdata[i+7*Q],
                                   qdata[i+8*Q],
                                   qdata[i+9*Q]
                                 };
    // -- Grad-to-Grad qdata
    const CeedScalar wBBJ[6]  =  { qdata[i+10*Q],
                                   qdata[i+11*Q],
                                   qdata[i+12*Q],
                                   qdata[i+13*Q],
                                   qdata[i+14*Q],
                                   qdata[i+15*Q]
                                 };  

    //--------- dU/dX
    const CeedScalar drhodX[3]  = { drho[0] * wBJ[0] + drho[1] * wBJ[1] + drho[2] * wBJ[2],
                                    drho[0] * wBJ[3] + drho[1] * wBJ[4] + drho[2] * wBJ[5],
                                    drho[0] * wBJ[6] + drho[1] * wBJ[7] + drho[2] * wBJ[8]
                                 };   
    const CeedScalar dUdX[3][3] ={ {dU[0][0] * wBJ[0] + dU[0][1] * wBJ[1] + dU[0][2] * wBJ[2],
                                    dU[0][0] * wBJ[3] + dU[0][1] * wBJ[4] + dU[0][2] * wBJ[5],
                                    dU[0][0] * wBJ[6] + dU[0][1] * wBJ[7] + dU[0][2] * wBJ[8]},
                                   {dU[1][0] * wBJ[0] + dU[1][1] * wBJ[1] + dU[1][2] * wBJ[2],
                                    dU[1][0] * wBJ[3] + dU[1][1] * wBJ[4] + dU[1][2] * wBJ[5],
                                    dU[1][0] * wBJ[6] + dU[1][1] * wBJ[7] + dU[1][2] * wBJ[8]},
                                   {dU[2][0] * wBJ[0] + dU[2][1] * wBJ[1] + dU[2][2] * wBJ[2], 
                                    dU[2][0] * wBJ[3] + dU[2][1] * wBJ[4] + dU[2][2] * wBJ[5],
                                    dU[2][0] * wBJ[6] + dU[2][1] * wBJ[7] + dU[2][2] * wBJ[8]}
                                  };
    const CeedScalar dEdX[3]  = { dE[0] * wBJ[0] + dE[1] * wBJ[1] + dE[2] * wBJ[2],
                                  dE[0] * wBJ[3] + dE[1] * wBJ[4] + dE[2] * wBJ[5],
                                  dE[0] * wBJ[6] + dE[1] * wBJ[7] + dE[2] * wBJ[8]
                                 }; 

    // ke = kinetic energy
    const CeedScalar ke = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/2;
    // P = pressure
    const CeedScalar P  =  ( E - ke * rho ) * (gamma - 1);
    // Tau = [TauC, TauM, TauM, TauM, TauE]
    const CeedScalar dt = 0.0001;        //1.e-5;
    const CeedScalar uiujgij = ( wBBJ[0]*u[0]*u[0] + wBBJ[3]*u[1]*u[1] + wBBJ[5]*u[2]*u[2] + 
                                 2*wBBJ[1]*u[0]*u[1] + 2*wBBJ[2]*u[0]*u[2] + 2*wBBJ[4]*u[1]*u[2])/wJ;
    const CeedScalar C1 =1.;
    const CeedScalar Cc =1.;
    const CeedScalar Ce =1.; 
    const CeedScalar f1   = rho * sqrt( 2/(C1 * dt) + uiujgij );
    const CeedScalar TauC = (Cc * f1 * wJ) / (8*(wBBJ[0] + wBBJ[3] + wBBJ[5]));      
    const CeedScalar TauM = 1/f1;
    const CeedScalar TauE = TauM/(Ce * cv);
    
    //-- Stabilizing terms
    // --- Stab[5][3] = Ai^T * Tau * Aj * U,j
    const CeedScalar S_conv[5][3] ={{TauM*(u[0]*u[0] - (Rd*ke)/cv)*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] -
                                     dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] +
                                     drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) + 
                                     TauM*u[0]*u[1]*(drho[1]*(u[1]*u[1] - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] -
                                     dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv +
                                     (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv) + TauM*u[0]*u[2]*(drho[2]*(u[2]*u[2]  - 
                                     (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + 
                                     drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv) + 
                                     TauE*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[2]*u[2] + ke))/cv) - dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) -
                                     dU[0][0]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) +
                                     drho[2]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv +
                                     (Rd*dU[1][0]*u[0]*u[1])/cv + (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + 
                                     (Rd*dU[2][1]*u[1]*u[2])/cv),
                                  
                                     TauM*(u[1]*u[1] - (Rd*ke)/cv)*(drho[1]*(u[1]*u[1] - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] -
                                     dU[2][2]*u[1] - dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] -
                                     (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv) + TauM*u[0]*u[1]*(drho[0]*(u[0]*u[0] -
                                     (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) +
                                     drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) +
                                     TauM*u[1]*u[2]*(drho[2]*(u[2]*u[2]  - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] -
                                     dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + 
                                     (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv) + TauE*u[1]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - dU[1][1]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[1]*u[1] + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - dE[0]*u[0]*(Rd/cv + 1) -
                                     dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv) +
                                     drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + 
                                     (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv),
                                     
                                     TauM*(u[2]*u[2] - (Rd*ke)/cv)*(drho[2]*(u[2]*u[2] - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - 
                                     dU[2][1]*u[1] - dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - 
                                     (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv) + TauM*u[0]*u[2]*(drho[0]*(u[0]*u[0] -
                                     (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) +
                                     drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) + 
                                     TauM*u[1]*u[2]*(drho[1]*(u[1]*u[1]  - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] - 
                                     dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + 
                                     (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv) + TauE*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - dU[1][1]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[1]*u[1] + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv)}, 
                                   
                                    {TauC*(dU[0][0] + dU[1][1] + dU[2][2]) - TauM*u[1]*(drho[1]*(u[1]*u[1]  - (Rd*ke)/cv) - dU[1][0]*u[0] - 
                                     dU[1][2]*u[2] - dU[2][2]*u[1] - dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + 
                                     drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv) - 
                                     TauM*u[2]*(drho[2]*(u[2]*u[2]  - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + 
                                     dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + 
                                     (Rd*dU[1][2]*u[1])/cv) - TauE*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv)*(drho[0]*u[0]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) - dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[2]*u[2]  + ke))/cv) - dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - 
                                     dU[0][0]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + 
                                     drho[2]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + 
                                     (Rd*dU[1][0]*u[0]*u[1])/cv + (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + 
                                     (Rd*dU[2][1]*u[1]*u[2])/cv) + TauM*u[0]*(Rd/cv - 2)*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - 
                                     dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + 
                                     drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv),
                                    
                                     (Rd*TauM*u[0]*(drho[1]*(u[1]*u[1] - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] - dU[0][0]*u[1] + 
                                     dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + 
                                     (Rd*dU[2][1]*u[2])/cv))/cv - TauM*u[1]*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - 
                                     dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - 
                                     (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) + (Rd*TauE*u[0]*u[1]*(drho[0]*u[0]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) - dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[2]*u[2]  + ke))/cv) - dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - 
                                     dU[0][0]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + 
                                     drho[2]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + 
                                     (Rd*dU[1][0]*u[0]*u[1])/cv + (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + 
                                     (Rd*dU[2][1]*u[1]*u[2])/cv))/cv,
                                    
                                     (Rd*TauM*u[0]*(drho[2]*(u[2]*u[2] - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + 
                                     dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + 
                                     (Rd*dU[1][2]*u[1])/cv))/cv - TauM*u[2]*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - 
                                     dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + 
                                     (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) + (Rd*TauE*u[0]*u[2]*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv))/cv},
                                   
                                    {(Rd*TauM*u[1]*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + 
                                     dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + 
                                     (Rd*dU[2][0]*u[2])/cv))/cv - TauM*u[0]*(drho[1]*(u[1]*u[1]  - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - 
                                     dU[2][2]*u[1] - dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + 
                                     (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv) + (Rd*TauE*u[0]*u[1]*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv))/cv,
                                   
                                     TauC*(dU[0][0] + dU[1][1] + dU[2][2]) - TauM*u[0]*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - 
                                     dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - 
                                     (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) - TauM*u[2]*(drho[2]*(u[2]*u[2] - (Rd*ke)/cv) - 
                                     dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + 
                                     drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv) - TauE*(E*(Rd/cv + 1) - 
                                     (Rd*(u[1]*u[1]  + ke))/cv)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - dU[1][1]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - dE[0]*u[0]*(Rd/cv + 1) - 
                                     dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv) + 
                                     drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + 
                                     (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + (Rd*dU[1][2]*u[1]*u[2])/cv + 
                                     (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv) + TauM*u[1]*(Rd/cv - 2)*(drho[1]*(u[1]*u[1]  - (Rd*ke)/cv) - 
                                     dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] - dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + 
                                     drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv),
                                    
                                     (Rd*TauM*u[1]*(drho[2]*(u[2]*u[2] - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + 
                                     dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + 
                                     (Rd*dU[1][2]*u[1])/cv))/cv - TauM*u[2]*(drho[1]*(u[1]*u[1]  - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - 
                                     dU[2][2]*u[1] - dU[0][0]*u[1] + dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + 
                                     (Rd*dU[0][1]*u[0])/cv + (Rd*dU[2][1]*u[2])/cv) + (Rd*TauE*u[1]*u[2]*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv))/cv}, 
                                   
                                    {(Rd*TauM*u[2]*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + 
                                     dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + 
                                     (Rd*dU[2][0]*u[2])/cv))/cv - TauM*u[0]*(drho[2]*(u[2]*u[2]  - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - 
                                     dU[2][1]*u[1] - dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + 
                                     (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv) + (Rd*TauE*u[0]*u[2]*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv))/cv,
                                   
                                     (Rd*TauM*u[2]*(drho[1]*(u[1]*u[1] - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] - dU[0][0]*u[1] + 
                                     dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + 
                                     (Rd*dU[2][1]*u[2])/cv))/cv - TauM*u[1]*(drho[2]*(u[2]*u[2]  - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - 
                                     dU[2][1]*u[1] - dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + 
                                     (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv) + (Rd*TauE*u[1]*u[2]*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv))/cv,
                                   
                                     TauC*(dU[0][0] + dU[1][1] + dU[2][2]) - TauM*u[0]*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - 
                                     dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + 
                                     drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + (Rd*dU[2][0]*u[2])/cv) - 
                                     TauM*u[1]*(drho[1]*(u[1]*u[1]  - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] - dU[0][0]*u[1] + 
                                     dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + 
                                     (Rd*dU[2][1]*u[2])/cv) - TauE*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv)*(drho[0]*u[0]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) - dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[2]*u[2]  + ke))/cv) - dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - 
                                     dU[0][0]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + 
                                     drho[2]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + 
                                     (Rd*dU[1][0]*u[0]*u[1])/cv + (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + 
                                     (Rd*dU[2][1]*u[1]*u[2])/cv) + TauM*u[2]*(Rd/cv - 2)*(drho[2]*(u[2]*u[2]  - (Rd*ke)/cv) - dU[1][1]*u[2] - 
                                     dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + 
                                     drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + (Rd*dU[1][2]*u[1])/cv)}, 
                                   
                                    {-(Rd*TauM*(drho[0]*(u[0]*u[0] - (Rd*ke)/cv) - dU[0][2]*u[2] - dU[1][1]*u[0] - dU[2][2]*u[0] - dU[0][1]*u[1] + 
                                     dU[0][0]*u[0]*(Rd/cv - 2) + drho[1]*u[0]*u[1] + drho[2]*u[0]*u[2] - (Rd*dE[0])/cv + (Rd*dU[1][0]*u[1])/cv + 
                                     (Rd*dU[2][0]*u[2])/cv))/cv - TauE*u[0]*(Rd/cv + 1)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv),
                                   
                                     -(Rd*TauM*(drho[1]*(u[1]*u[1] - (Rd*ke)/cv) - dU[1][0]*u[0] - dU[1][2]*u[2] - dU[2][2]*u[1] - dU[0][0]*u[1] + 
                                     dU[1][1]*u[1]*(Rd/cv - 2) + drho[0]*u[0]*u[1] + drho[2]*u[1]*u[2] - (Rd*dE[1])/cv + (Rd*dU[0][1]*u[0])/cv + 
                                     (Rd*dU[2][1]*u[2])/cv))/cv - TauE*u[1]*(Rd/cv + 1)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv),
                                   
                                     -(Rd*TauM*(drho[2]*(u[2]*u[2] - (Rd*ke)/cv) - dU[1][1]*u[2] - dU[2][0]*u[0] - dU[2][1]*u[1] - dU[0][0]*u[2] + 
                                     dU[2][2]*u[2]*(Rd/cv - 2) + drho[0]*u[0]*u[2] + drho[1]*u[1]*u[2] - (Rd*dE[2])/cv + (Rd*dU[0][2]*u[0])/cv + 
                                     (Rd*dU[1][2]*u[1])/cv))/cv - TauE*u[2]*(Rd/cv + 1)*(drho[0]*u[0]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) - 
                                     dU[1][1]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv) - dU[2][2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - 
                                     dE[0]*u[0]*(Rd/cv + 1) - dE[1]*u[1]*(Rd/cv + 1) - dE[2]*u[2]*(Rd/cv + 1) - dU[0][0]*(E*(Rd/cv + 1) - 
                                     (Rd*(u[0]*u[0] + ke))/cv) + drho[1]*u[1]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv) + drho[2]*u[2]*(E*(Rd/cv + 1) - 
                                     (2*Rd*ke)/cv) + (Rd*dU[0][1]*u[0]*u[1])/cv + (Rd*dU[0][2]*u[0]*u[2])/cv + (Rd*dU[1][0]*u[0]*u[1])/cv + 
                                     (Rd*dU[1][2]*u[1]*u[2])/cv + (Rd*dU[2][0]*u[0]*u[2])/cv + (Rd*dU[2][1]*u[1]*u[2])/cv)}
                                  };                            
    //-----Body Force Ai^T * Tau * Fb
    const CeedScalar S_Fb[5][3] =  { { -TauM*g*rho*u[0]*u[2] - TauE*g*rho*u[0]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv),
                                       -TauM*g*rho*u[1]*u[2] - TauE*g*rho*u[1]*u[2]*(E*(Rd/cv + 1) - (2*Rd*ke)/cv),
                                       -TauM*g*rho*(u[2]*u[2]  - (Rd*ke)/cv) - TauE*g*rho*u[2]*u[2] *(E*(Rd/cv + 1) - (2*Rd*ke)/cv) }, 
                                     { TauM*g*rho*u[2] + TauE*g*rho*u[2]*(E*(Rd/cv + 1) - (Rd*(u[0]*u[0] + ke))/cv),
                                       -(Rd*TauE*g*rho*u[0]*u[1]*u[2])/cv, 
                                       -(Rd*TauM*g*rho*u[0])/cv - (Rd*TauE*g*rho*u[0]*u[2]*u[2] )/cv },
                                     { -(Rd*TauE*g*rho*u[0]*u[1]*u[2])/cv,
                                       TauM*g*rho*u[2] + TauE*g*rho*u[2]*(E*(Rd/cv + 1) - (Rd*(u[1]*u[1]  + ke))/cv),
                                       -(Rd*TauM*g*rho*u[1])/cv - (Rd*TauE*g*rho*u[1]*u[2]*u[2] )/cv }, 
                                     { TauM*g*rho*u[0] - (Rd*TauE*g*rho*u[0]*u[2]*u[2])/cv,
                                       TauM*g*rho*u[1] - (Rd*TauE*g*rho*u[1]*u[2]*u[2])/cv,
                                       TauE*g*rho*u[2]*(E*(Rd/cv + 1) - (Rd*(u[2]*u[2]  + ke))/cv) - TauM*g*rho*u[2]*(Rd/cv - 2) }, 
                                     { TauE*g*rho*u[0]*u[2]*(Rd/cv + 1),
                                       TauE*g*rho*u[1]*u[2]*(Rd/cv + 1),
                                       TauE*g*rho*(Rd/cv + 1)*u[2]*u[2]  + (Rd*TauM*g*rho)/cv } 
                                   };
    // The Physics
    // -- Density--- u rho
    dv[i+(0+5*0)*Q]  = rho*u[0]*wBJ[0] + rho*u[1]*wBJ[1] + rho*u[2]*wBJ[2];
    dv[i+(0+5*1)*Q]  = rho*u[0]*wBJ[3] + rho*u[1]*wBJ[4] + rho*u[2]*wBJ[5];
    dv[i+(0+5*2)*Q]  = rho*u[0]*wBJ[6] + rho*u[1]*wBJ[7] + rho*u[2]*wBJ[8];
    // -- Momentum-- rho (u x u) + P I3
    dv[i+(1+5*0)*Q]  = (P + rho*u[0]*u[0])*wBJ[0] + rho*u[0]*u[1] *wBJ[1] + rho*u[0]*u[2] *wBJ[2];
    dv[i+(1+5*1)*Q]  = (P + rho*u[0]*u[0])*wBJ[3] + rho*u[0]*u[1] *wBJ[4] + rho*u[0]*u[2]* wBJ[5];
    dv[i+(1+5*2)*Q]  = (P + rho*u[0]*u[0])*wBJ[6] + rho*u[0]*u[1] *wBJ[7] + rho*u[0]*u[2]* wBJ[8];
    dv[i+(2+5*0)*Q]  = rho*u[1]*u[0]* wBJ[0] + (P + rho*u[1]*u[1])*wBJ[1] + rho*u[1]*u[2]* wBJ[2];
    dv[i+(2+5*1)*Q]  = rho*u[1]*u[0]* wBJ[3] + (P + rho*u[1]*u[1])*wBJ[4] + rho*u[1]*u[2]* wBJ[5];  
    dv[i+(2+5*2)*Q]  = rho*u[1]*u[0]* wBJ[6] + (P + rho*u[1]*u[1])*wBJ[7] + rho*u[1]*u[2]* wBJ[8];
    dv[i+(3+5*0)*Q]  = rho*u[2]*u[0]* wBJ[0] + rho*u[2]*u[1]* wBJ[1] + (P +rho*u[2]*u[2])* wBJ[2];
    dv[i+(3+5*1)*Q]  = rho*u[2]*u[0]* wBJ[3] + rho*u[2]*u[1]* wBJ[4] + (P +rho*u[2]*u[2])* wBJ[5];
    dv[i+(3+5*2)*Q]  = rho*u[2]*u[0]* wBJ[6] + rho*u[2]*u[1]* wBJ[7] + (P +rho*u[2]*u[2])* wBJ[8];
    // -- Total Energy-- (E + P) u
    dv[i+(4+5*0)*Q]  = (E + P)*(u[0]*wBJ[0] + u[1]*wBJ[1] + u[2]*wBJ[2]);
    dv[i+(4+5*1)*Q]  = (E + P)*(u[0]*wBJ[3] + u[1]*wBJ[4] + u[2]*wBJ[5]);
    dv[i+(4+5*2)*Q]  = (E + P)*(u[0]*wBJ[6] + u[1]*wBJ[7] + u[2]*wBJ[8]);
    //----------Body Force (F)
    v[i+0*Q] = 0;
    v[i+1*Q] = 0;
    v[i+2*Q] = 0;
    v[i+3*Q] = rho*g      * wJ;
    v[i+4*Q] = rho*g*u[2] * wJ;
    // ---- Convective Stabilizing Terms
    dv[i+(1+5*0)*Q] -= S_conv[0][0] * wBJ[0] + S_conv[0][1] * wBJ[1] + S_conv[0][2] * wBJ[2];
    dv[i+(1+5*1)*Q] -= S_conv[0][0] * wBJ[3] + S_conv[0][1] * wBJ[4] + S_conv[0][2] * wBJ[5];
    dv[i+(1+5*2)*Q] -= S_conv[0][0] * wBJ[6] + S_conv[0][1] * wBJ[7] + S_conv[0][2] * wBJ[8];

    dv[i+(1+5*0)*Q] -= S_conv[1][0] * wBJ[0] + S_conv[1][1] * wBJ[1] + S_conv[1][2] * wBJ[2];
    dv[i+(1+5*1)*Q] -= S_conv[1][0] * wBJ[3] + S_conv[1][1] * wBJ[4] + S_conv[1][2] * wBJ[5];
    dv[i+(1+5*2)*Q] -= S_conv[1][0] * wBJ[6] + S_conv[1][1] * wBJ[7] + S_conv[1][2] * wBJ[8];

    dv[i+(2+5*0)*Q] -= S_conv[2][0] * wBJ[0] + S_conv[2][1] * wBJ[1] + S_conv[2][2] * wBJ[2];
    dv[i+(2+5*1)*Q] -= S_conv[2][0] * wBJ[3] + S_conv[2][1] * wBJ[4] + S_conv[2][2] * wBJ[5];
    dv[i+(2+5*2)*Q] -= S_conv[2][0] * wBJ[6] + S_conv[2][1] * wBJ[7] + S_conv[2][2] * wBJ[8];

    dv[i+(3+5*0)*Q] -= S_conv[3][0] * wBJ[0] + S_conv[3][1] * wBJ[1] + S_conv[3][2] * wBJ[2];
    dv[i+(3+5*1)*Q] -= S_conv[3][0] * wBJ[3] + S_conv[3][1] * wBJ[4] + S_conv[3][2] * wBJ[5];
    dv[i+(3+5*2)*Q] -= S_conv[3][0] * wBJ[6] + S_conv[3][1] * wBJ[7] + S_conv[3][2] * wBJ[8];

    dv[i+(4+5*0)*Q] -= S_conv[4][0] * wBJ[0] + S_conv[4][1] * wBJ[1] + S_conv[4][2] * wBJ[2];
    dv[i+(4+5*1)*Q] -= S_conv[4][0] * wBJ[3] + S_conv[4][1] * wBJ[4] + S_conv[4][2] * wBJ[5];
    dv[i+(4+5*2)*Q] -= S_conv[4][0] * wBJ[6] + S_conv[4][1] * wBJ[7] + S_conv[4][2] * wBJ[8];
    // ---- Body Force Stabilizing Terms
    dv[i+(1+5*0)*Q] += S_Fb[0][0] * wBJ[0] + S_Fb[0][1] * wBJ[1] + S_Fb[0][2] * wBJ[2];
    dv[i+(1+5*1)*Q] += S_Fb[0][0] * wBJ[3] + S_Fb[0][1] * wBJ[4] + S_Fb[0][2] * wBJ[5];
    dv[i+(1+5*2)*Q] += S_Fb[0][0] * wBJ[6] + S_Fb[0][1] * wBJ[7] + S_Fb[0][2] * wBJ[8];

    dv[i+(1+5*0)*Q] += S_Fb[1][0] * wBJ[0] + S_Fb[1][1] * wBJ[1] + S_Fb[1][2] * wBJ[2];
    dv[i+(1+5*1)*Q] += S_Fb[1][0] * wBJ[3] + S_Fb[1][1] * wBJ[4] + S_Fb[1][2] * wBJ[5];
    dv[i+(1+5*2)*Q] += S_Fb[1][0] * wBJ[6] + S_Fb[1][1] * wBJ[7] + S_Fb[1][2] * wBJ[8];

    dv[i+(2+5*0)*Q] += S_Fb[2][0] * wBJ[0] + S_Fb[2][1] * wBJ[1] + S_Fb[2][2] * wBJ[2];
    dv[i+(2+5*1)*Q] += S_Fb[2][0] * wBJ[3] + S_Fb[2][1] * wBJ[4] + S_Fb[2][2] * wBJ[5];
    dv[i+(2+5*2)*Q] += S_Fb[2][0] * wBJ[6] + S_Fb[2][1] * wBJ[7] + S_Fb[2][2] * wBJ[8];

    dv[i+(3+5*0)*Q] += S_Fb[3][0] * wBJ[0] + S_Fb[3][1] * wBJ[1] + S_Fb[3][2] * wBJ[2];
    dv[i+(3+5*1)*Q] += S_Fb[3][0] * wBJ[3] + S_Fb[3][1] * wBJ[4] + S_Fb[3][2] * wBJ[5];
    dv[i+(3+5*2)*Q] += S_Fb[3][0] * wBJ[6] + S_Fb[3][1] * wBJ[7] + S_Fb[3][2] * wBJ[8];

    dv[i+(4+5*0)*Q] += S_Fb[4][0] * wBJ[0] + S_Fb[4][1] * wBJ[1] + S_Fb[4][2] * wBJ[2];
    dv[i+(4+5*1)*Q] += S_Fb[4][0] * wBJ[3] + S_Fb[4][1] * wBJ[4] + S_Fb[4][2] * wBJ[5];
    dv[i+(4+5*2)*Q] += S_Fb[4][0] * wBJ[6] + S_Fb[4][1] * wBJ[7] + S_Fb[4][2] * wBJ[8];    

  } // End Quadrature Point Loop

  // Return
  return 0;
}

// *****************************************************************************
#endif
