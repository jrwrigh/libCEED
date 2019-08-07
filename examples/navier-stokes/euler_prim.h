/// @file
/// Density current initial condition and operator for Navier-Stokes example using PETSc

// Model from:
//   Semi-Implicit Formulations of the Navier-Stokes Equations: Application to
//   Nonhydrostatic Atmospheric Modeling, Giraldo, Restelli, and Lauter (2010).

#ifndef euler_prim_h
#define euler_prim_h

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
//   Exner pressure. Initial velocity is zero.
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
//   
//
//  Boundary Conditions:
//       Velocity:
//              ui = 0
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
// *****************************************************************************
static int ICsDCprim(void *ctx, CeedInt Q,
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


    // Initial Conditions
    q0[i+0*Q] = Pi;
    q0[i+1*Q] = 0.0;
    q0[i+2*Q] = 0.0;
    q0[i+3*Q] = 0.0;
    q0[i+4*Q] = theta;

    // Homogeneous Dirichlet Boundary Conditions for Velocity
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
// This QFunction implements the following formulation of Euler Equation.
//
// This is Euler Equation in primitive form with state
//   variables of pressure, velocity, and temperature.
//
// Euler Equation:
//   drho/dt + div( U )                               = 0
//   dU/dt   + div( rho (u x u) + P I3 ) + rho g khat = 0
//   dE/dt   + div( (E + P) u )                       = 0
//
//----Changing to primitive variables:
//    Y,t + {A0}^-1 * Ai * Y,i - Body Force =0
//
//     State Variables: q = ( P, u1, u2, u3, T )
//
// Constants:
//   cv              ,  Specific heat, constant volume
//   cp              ,  Specific heat, constant pressure
//   g               ,  Gravity
//   Rd     = cp - cv,  Specific heat difference
// *****************************************************************************
static int DCprim(void *ctx, CeedInt Q,
              const CeedScalar *const *in, CeedScalar *const *out) {
  // Inputs
  const CeedScalar *q = in[0], *dq = in[1], *qdata = in[2], *x = in[3];
  // Outputs
  CeedScalar *v = out[0], *dv = out[1];
  // Context
  const CeedScalar *context = (const CeedScalar*)ctx;
  const CeedScalar cv         = context[3];
  const CeedScalar g          = context[5];
  const CeedScalar Rd         = context[6];

  CeedPragmaOMP(simd)
  // Quadrature Point Loop
  for (CeedInt i=0; i<Q; i++) {

    //--------------------------------SETUP---------------------------------//
    // --------- Interp in
    const CeedScalar P       =   q[i+0*Q];
    
    const CeedScalar u[3]    = { q[i+1*Q],
                                 q[i+2*Q],
                                 q[i+3*Q]
                               };
    
    const CeedScalar T       = q[i+4*Q];  

    // -- rho
    const CeedScalar rho     =   P/(Rd*T);

    // --------- Grad in
    const CeedScalar dP[3]   =  { dq[i+(0+5*0)*Q],
                                  dq[i+(0+5*1)*Q],
                                  dq[i+(0+5*2)*Q]
                             };
    
    const CeedScalar du[3][3]= { {dq[i+(1+5*0)*Q], 
                                  dq[i+(1+5*1)*Q], 
                                  dq[i+(1+5*2)*Q]},  
                                 {dq[i+(2+5*0)*Q],
                                  dq[i+(2+5*1)*Q], 
                                  dq[i+(2+5*2)*Q]}, 
                                 {dq[i+(3+5*0)*Q], 
                                  dq[i+(3+5*1)*Q], 
                                  dq[i+(3+5*2)*Q]}  
                               };
 
    const CeedScalar dT[3]   = {  dq[i+(4+5*0)*Q],
                                  dq[i+(4+5*1)*Q],
                                  dq[i+(4+5*2)*Q]
                               };



    // --- Interp-to-Interp qdata
    const CeedScalar wJ       =   qdata[i+ 0*Q];     


    // --- Interp-to-Grad qdata
    const CeedScalar wBJ[9]   = { qdata[i+ 1*Q],
                                  qdata[i+ 2*Q],
                                  qdata[i+ 3*Q],
                                  qdata[i+ 4*Q],
                                  qdata[i+ 5*Q],
                                  qdata[i+ 6*Q],
                                  qdata[i+ 7*Q],
                                  qdata[i+ 8*Q],
                                  qdata[i+ 9*Q]
                                };

    //------------ dY/dX
    const CeedScalar dP_dX[3]  = { dP[0] * wBJ[0] + dP[1] * wBJ[1] + dP[2] * wBJ[2],
                                   dP[0] * wBJ[3] + dP[1] * wBJ[4] + dP[2] * wBJ[5],
                                   dP[0] * wBJ[6] + dP[1] * wBJ[7] + dP[2] * wBJ[8]
                                 };


    const CeedScalar du_dX[3][3] ={ {du[0][0] * wBJ[0] + du[0][1] * wBJ[1] + du[0][2] * wBJ[2],
                                     du[0][0] * wBJ[3] + du[0][1] * wBJ[4] + du[0][2] * wBJ[5],
                                     du[0][0] * wBJ[6] + du[0][1] * wBJ[7] + du[0][2] * wBJ[8]},
                                    {du[1][0] * wBJ[0] + du[1][1] * wBJ[1] + du[1][2] * wBJ[2],
                                     du[1][0] * wBJ[3] + du[1][1] * wBJ[4] + du[1][2] * wBJ[5],
                                     du[1][0] * wBJ[6] + du[1][1] * wBJ[7] + du[1][2] * wBJ[8]},
                                    {du[2][0] * wBJ[0] + du[2][1] * wBJ[1] + du[2][2] * wBJ[2], 
                                     du[2][0] * wBJ[3] + du[2][1] * wBJ[4] + du[2][2] * wBJ[5],
                                     du[2][0] * wBJ[6] + du[2][1] * wBJ[7] + du[2][2] * wBJ[8]}
                                  };


    const CeedScalar dT_dX[3]  = { dT[0] * wBJ[0] + dT[1] * wBJ[1] + dT[2] * wBJ[2],
                                   dT[0] * wBJ[3] + dT[1] * wBJ[4] + dT[2] * wBJ[5],
                                   dT[0] * wBJ[6] + dT[1] * wBJ[7] + dT[2] * wBJ[8]
                                 };                    
                                                 
    //------------------------------>  THE PHYSICS  <------------------------------//

    // --------Convective Flux
              //A0inv * Ai * (Y,i * wBJ)                  

    dv[i+(0+5*0)*Q]  = 0;
    dv[i+(0+5*1)*Q]  = 0;
    dv[i+(0+5*2)*Q]  = 0;

    dv[i+(1+5*0)*Q]  = 0;
    dv[i+(1+5*1)*Q]  = 0;
    dv[i+(1+5*2)*Q]  = 0;

    dv[i+(2+5*0)*Q]  = 0;
    dv[i+(2+5*1)*Q]  = 0;
    dv[i+(2+5*2)*Q]  = 0;

    dv[i+(3+5*0)*Q]  = 0;
    dv[i+(3+5*1)*Q]  = 0; 
    dv[i+(3+5*2)*Q]  = 0;

    dv[i+(4+5*0)*Q]  = 0;
    dv[i+(4+5*1)*Q]  = 0;
    dv[i+(4+5*2)*Q]  = 0;
                      

    //---------Non-divergence terms
             // A0inv * (Ai * Y,i - F)

    v[i+0*Q] = u[0] * dP_dX[0] + (Rd*Rd*T*rho*(cv/Rd + 1)/cv) * du_dX[0][0] +
               u[1] * dP_dX[1] + (Rd*Rd*T*rho*(cv/Rd + 1)/cv) * du_dX[1][1] +
               u[2] * dP_dX[2] + (Rd*Rd*T*rho*(cv/Rd + 1)/cv) * du_dX[2][2]; 
               
    v[i+1*Q] = u[0] * du_dX[0][0] + dP_dX[0]/rho +
               u[1] * du_dX[0][1] +
               u[2] * du_dX[0][2];

    v[i+2*Q] = u[0] * du_dX[1][0] +
               u[1] * du_dX[1][1] + dP_dX[1]/rho +
               u[2] * du_dX[1][2];

    v[i+3*Q] = u[0] * du_dX[2][0] +
               u[1] * du_dX[2][1] +
               u[2] * du_dX[2][2] + dP_dX[2]/rho -
               g * wJ; 

    v[i+4*Q] = u[0] * dT_dX[0] + (Rd*T/cv) * du_dX[0][0] +
               u[1] * dT_dX[1] + (Rd*T/cv) * du_dX[1][1] +
               u[2] * dT_dX[2] + (Rd*T/cv) * du_dX[2][2];      


  } // End Quadrature Point Loop

  // Return
  return 0;
}

// *****************************************************************************
#endif
