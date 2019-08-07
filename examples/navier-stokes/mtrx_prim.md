````
These matrices are for primitive variables
*********************************************************************

ke = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/2               //kinetic energy 
---------------------------------------

A0[5][5] = {0}

A0[0][0] = 1/(Rd*T)
A0[0][4] = -rho/T

A0[1][0] = u[0]/(Rd*T)
A0[1][2] = rho
A0[1][4] = -rho*u[0]/T

A0[2][0] = u[1]/(Rd*T)
A0[2][2] = rho
A0[2][4] = -rho*u[1]/T

A0[3][0] = u[2]/(Rd*T)
A0[3][3] = rho
A0[3][4] = -rho*u[2]/T

A0[4][0] = (cv + ke/T)/Rd
A0[4][1] = rho*u[0]
A0[4][2] = rho*u[1]
A0[4][3] = rho*u[2]
A0[4][4] = -rho*ke/T

****************

A0_inv[5][5] ={0}

A0_inv[0][0] = Rd*ke/cv
A0_inv[0][1] = -Rd*u[0]/cv
A0_inv[0][2] = -Rd*u[1]/cv
A0_inv[0][3] = -Rd*u[2]/cv
A0_inv[0][4] = Rd/cv

A0_inv[1][0] = -u[0]/rho
A0_inv[1][1] = 1/rho

A0_inv[2][0] = -u[1]/rho
A0_inv[2][2] = 1/rho

A0_inv[3][0] = -u[2]/rho
A0_inv[3][3] = 1/rho

A0_inv[4][0] = (ke - cv*T)/(cv*rho)
A0_inv[4][1] = -u[0]/(cv*rho)
A0_inv[4][2] = -u[1]/(cv*rho)
A0_inv[4][3] = -u[2]/(cv*rho)
A0_inv[4][4] =  1/(cv*rho)

*********************************************************************

Ai1bar = 1 + cv/Rd + ke/(Rd*T)
Ai2bar = rho*Rd*T*(1 + cv/Rd) + rho*ke
--------------------------------------------

Ai_1[5][5] = {0}

Ai_1[0][0] = u[0]/(Rd*T)
Ai_1[0][1] = rho
Ai_1[0][4] = -rho*u[0]/T

Ai_1[1][0] = u[0]*u[0]/(Rd*T) + 1
Ai_1[1][1] = 2*rho*u[0]
Ai_1[1][4] = -rho*u[0]*u[0]/T

Ai_1[2][1] = u[0]*u[1]/(Rd*T)
Ai_1[2][2] = rho*u[1]
Ai_1[2][3] = rho*u[0]
Ai_1[2][5] = -rho*u[0]*u[1]/T

Ai_1[3][0] = u[0]*u[2]/(Rd*T)
Ai_1[3][1] = rho*u[2]
Ai_1[3][3] = rho*u[0]
Ai_1[3][4] = -rho*u[0]*u[2]/T

Ai_1[4][0] = u[0]*Ai1bar
Ai_1[4][1] = Ai2bar + rho*u[0]*u[0]
Ai_1[4][2] = rho*u[0]*u[1]
Ai_1[4][3] = rho*u[0]*u[2]
Ai_1[4][4] = -rho*u[0]*ke/T

****************

Ai_2[5][5] = {0}

Ai_2[0][0] = u[1]/(Rd*T)
Ai_2[0][2] = rho
Ai_2[0][4] = -rho*u[1]/T

Ai_2[1][0] = u[1]*u[0]/(Rd*T)
Ai_2[1][1] = rho*u[1]
Ai_2[1][2] = rho*u[0]
Ai_2[1][4] = -rho*u[1]*u[0]/T

Ai_2[2][0] = u[1]*u[1]/(Rd*T) +1
Ai_2[2][2] = 2*rho*u[1]
Ai_2[2][4] = -rho*u[1]*u[1]/T

Ai_2[3][0] = u[1]*u[2]/(Rd*T)
Ai_2[3][2] = rho*u[2]
Ai_2[3][3] = rho*u[1]
Ai_2[3][4] = -rho*u[1]*u[2]/T

Ai_2[4][0] = u[1]*Ai1bar
Ai_2[4][1] = rho*u[1]*u[0]
Ai_2[4][2] = Ai2bar + rho*u[1]*u[1]
Ai_2[4][3] = rho*u[1]*u[2]
Ai_2[4][4] = -rho*u[1]*ke/T

****************

Ai_3[5][5] = {0}

Ai_3[0][0] = u[2]/(Rd*T)
Ai_3[0][3] = rho
Ai_3[0][4] = -rho*u[2]/T

Ai_3[1][0] = u[2]*u[0]/(Rd*T)
Ai_3[1][1] = rho*u[2]
Ai_3[1][3] = rho*u[0]
Ai_3[1][4] = -rho*u[2]*u[0]/T

Ai_3[2][0] = u[2]*u[1]/(Rd*T)
Ai_3[2][2] = rho*u[2]
Ai_3[2][3] = rho*u[1]
Ai_3[2][4] = -rho*u[2]*u[1]/T

Ai_3[3][0] = u[2]*u[2]/(Rd*T) +1
Ai_3[3][3] = 2*rho*u[2]
Ai_3[3][4] = -rho*u[2]*u[2]/T

Ai_3[4][0] = u[2]*Ai1bar
Ai_3[4][1] = rho*u[2]*u[0]
Ai_3[4][2] = rho*u[2]*u[1]
Ai_3[4][3] = Ai2bar + rho*u[2]*u[2]
Ai_3[4][4] = -rho*u[2]*ke/T
*********************************************************************

K_ij[15][15] = {0}

//i=1, j=1
K_ij[1][1]   = 2*mu + lambda
K_ij[2][2]   = mu
K_ij[3][3]   = mu
K_ij[4][1]   = (2*mu + lambda)*u[0]
K_ij[4][2]   = mu*u[1]
K_ij[4][3]   = mu*u[2]
K_ij[4][4]   = k


//i=1, j=2
K_ij[1][7]   = lambda
K_ij[2][6]   = mu
K_ij[4][6]   = mu*u[1]
K_ij[4][7]   = mu*u[0]


//i=1, j=3
K_ij[1][13]  = lambda
K_ij[3][11]  = mu
K_ij[4][11]  = mu*u[2]
K_ij[4][13]  = mu*u[0]


//i=2, j=1
K_ij[6][2]   = mu
K_ij[7][1]   = lambda
K_ij[9][1]  = lambda*u[1]
K_ij[9][2]  = mu*u[0]


//i=2, j=2
K_ij[6][6]   = mu
K_ij[7][7]   = 2*mu + lambda
K_ij[8][8]   = mu
K_ij[9][6]  = mu*u[0]
K_ij[9][7]  = (2*mu + lambda)*u[1]
K_ij[9][8]  = mu*u[2]
K_ij[9][9] = k


//i=2, j=3
K_ij[7][13]  = lambda
K_ij[8][12]  = mu
K_ij[9][12] = mu*u[2]
K_ij[9][13] = lambda*u[1]


//i=3, j=1
K_ij[11][3]  = mu
K_ij[13][1]  = lambda
K_ij[14][1]  = lambda*u[2]
K_ij[14][3]  = mu*u[0]


//i=3, j=2
K_ij[12][8]  = mu
K_ij[13][7]  = lambda
K_ij[14][7]  = lambda*u[2]
K_ij[14][8]  = mu*u[1]


//i=3, j=3
K_ij[11][11] = mu
K_ij[12][12] = mu
K_ij[13][13] = 2*mu + lambda
K_ij[14][11] = mu*u[0]
K_ij[14][12] = mu*u[1]
K_ij[14][13] = (2*mu + lambda)*u[2]
K_ij[14][14] = k

*********************************************************************

A0inv_A1[5][5] = {0}

A0inv_A1[0][0] = u[0]
A0inv_A1[0][1] = Rd**2*T*rho*(cv/Rd + 1)/cv

A0inv_A1[1][0] = 1/rho
A0inv_A1[1][1] = u[0]

A0inv_A1[2][2] = u[0]

A0inv_A1[3][3] = u[0]

A0inv_A1[4][1] = Rd*T/cv
A0inv_A1[4][4] = u[0]
****************

A0inv_A2[5][5] = {0}

A0inv_A2[0][0] = u[1]
A0inv_A2[0][2] = Rd**2*T*rho*(cv/Rd + 1)/cv
 
A0inv_A2[1][1] = u[1]

A0inv_A2[2][0] = 1/rho
A0inv_A2[2][2] = u[1]

A0inv_A2[3][3] = u[1]

A0inv_A2[4][2] = Rd*T/cv
A0inv_A2[4][4] = u[1]
****************
A0inv_A3[5][5] = {0}

A0inv_A3[0][0] = u[2]
A0inv_A3[0][3] = Rd**2*T*rho*(cv/Rd + 1)/cv
 
A0inv_A3[1][1] = u[2]

A0inv_A3[2][2] = u[2]

A0inv_A3[3][0] = 1/rho
A0inv_A3[3][3] = u[2]

A0inv_A3[4][3] = Rd*T/cv
A0inv_A3[4][4] = u[2]
*********************************************************************

//-------> A0inv * A1 * dY_dx
 [0] = u[0] * dP_dX[0]    + (Rd**2*T*rho*(cv/Rd + 1)/cv) * du_dX[0][0]
 [1] = u[0] * du_dX[0][0] + dP_dX[0]/rho
 [2] = u[0] * du_dX[1][0]
 [3] = u[0] * du_dX[2][0]
 [4] = u[0] * dT_dX[0]    + (Rd*T/cv) * du_dX[0][0]

//-------> A0inv * A2 * dY_dy

 [5] = u[1] * dP_dX[1]    + (Rd**2*T*rho*(cv/Rd + 1)/cv) * du_dX[1][1]
 [6] = u[1] * du_dX[0][1]
 [7] = u[1] * du_dX[1][1] + dP_dX[1]/rho
 [8] = u[1] * du_dX[2][1]
 [9] = u[1] * dT_dX[1]    + (Rd*T/cv) * du_dX[1][1]

//-------> A0inv * A3 * dY_dz
 [10] = u[2] * dP_dX[2]    + (Rd**2*T*rho*(cv/Rd + 1)/cv) * du_dX[2][2]
 [11] = u[2] * du_dX[0][2]
 [12] = u[2] * du_dX[1][2]
 [13] = u[2] * du_dX[2][2] + dP_dX[2]/rho
 [14] = u[2] * dT_dX[2]    + (Rd*T/cv) * du_dX[2][2]

````











