````
uu = u[1]*u[1] + u[2]*u[2] + u[3]*u[3]


A0[5][5] = {0}

A0[1][1] = 1/(Rd*T)
A0[1][5] = -rho/(T)

A0[2][1] = u[1]/(Rd*T)
A0[2][3] = rho
A0[2][5] = -rho*u[1]/T

A0[3][1] = u[2]/(Rd*T)
A0[3][3] = rho
A0[3][5] = -rho*u[2]/T

A0[4][1] = u[3]/(Rd*T)
A0[4][4] = rho
A0[4][5] = -rho*u[3]/T

A0[5][1] = (cv + uu/(2*T) + g*z/T - (cp*theta0*exp(N*N*z/g))/T)/Rd
A0[5][2] = rho*u[1]
A0[5][3] = rho*u[2]
A0[5][4] = rho*u[3]
A0[5][5] = rho*(-g*z/T - uu/(2*T) + g*gcexp(-N*N*z/g)/(theta0*N*N))
*********************************************************************

Aia = (g*z)/(Rd*T) - rho*cp*theta0
Aib = -rho*g*z/T + rho*g*g*exp(-N*N*z/g)/(theta0*N*N)
Aic = 1 + cv/Rd + g*z/Rd + uu/(2*Rd*T) - rho*cp*theta0*exp(N*N*z/g)
Aid = Pi*(1 + cv/Rd) + rho*uu/2 + rho*g*z
Aie = -g*z/T - uu/(2*T) + g*g*exp(-N*N*z/g)/(theta0*N*N)


Ai_1[5][5] = {0}

Ai_1[1][1] = u[1]/(Rd*T)
Ai_1[1][2] = rho
Ai_1[1][5] = -rho*u[1]/T

Ai_1[2][1] = Aia + u[1]*u[1]/(Rd*T) + 1
Ai_1[2][2] = 2*rho*u[1]
Ai_1[2][5] = Aib - rho*u[1]*u[1]/T

Ai_1[3][1] = Aia + u[1]*u[2]/(Rd*T)
Ai_1[3][2] = rho*u[2]
Ai_1[3][3] = rho*u[1]
Ai_1[2][5] = Aib - rho*u[1]*u[2]/T

Ai_1[4][1] = Aia + u[1]*u[3]/(Rd*T)
Ai_1[4][2] = rho*u[3]
Ai_1[4][4] = rho*u[1]
Ai_1[4][5] = Aib - rho*u[1]*u[3]/T

Ai_1[5][1] = u[1]*Aic
Ai_1[5][2] = Aid + rho*u[1]*u[1]
Ai_1[5][3] = rho*u[1]*u[2]
Ai_1[5][4] = rho*u[1]*u[3]
Ai_1[5][5] = rho*u[1]*Aie

****************

Ai_2[5][5] = {0}

Ai_2[1][1] = u[2]/(Rd*T)
Ai_2[1][3] = rho
Ai_2[1][5] = -rho*u[2]/T

Ai_2[2][1] = Aia + u[2]*u[1]/(Rd*T)
Ai_2[2][2] = rho*u[2]
Ai_2[2][3] = rho*u[1]
Ai_2[2][5] = Aib - rho*u[2]*u[1]/T

Ai_2[3][1] = Aia + u[2]*u[2]/(Rd*T) +1
Ai_2[3][3] = 2*rho*u[2]
Ai_2[3][5] = Aib - rho*u[2]*u[2]/T

Ai_2[4][1] = Aia + u[2]*u[3]/(Rd*T)
Ai_2[4][3] = rho*u[3]
Ai_2[4][4] = rho*u[2]
Ai_2[4][5] = Aib - rho*u[2]*u[3]/T

Ai_2[5][1] = u[2]*Aic
Ai_2[5][2] = rho*u[2]*u[1]
Ai_2[5][3] = Aid + rho*u[2]*u[2]
Ai_2[5][4] = rho*u[2]*u[3]
Ai_2[5][5] = rho*u[2]*Aie

****************

Ai_3[5][5] = {0}

Ai_3[1][1] = u[3]/(Rd*T)
Ai_3[1][4] = rho
Ai_3[1][5] = -rho*u[3]/T

Ai_3[2][1] = Aia + u[3]*u[1]/(Rd*T)
Ai_3[2][2] = rho*u[3]
Ai_3[2][4] = rho*u[1]
Ai_3[2][5] = Aib - rho*u[3]*u[1]/T

Ai_3[3][1] = Aia + u[3]*u[2]/(Rd*T)
Ai_3[3][3] = rho*u[3]
Ai_3[3][4] = rho*u[2]
Ai_3[3][5] = Aib - rho*u[3]*u[2]/T

Ai_3[4][1] = Aia + u[3]*u[3]/(Rd*T) +1
Ai_3[4][4] = 2*rho*u[3]
Ai_3[4][5] = Aib - rho*u[3]*u[3]/T

Ai_3[5][1] = u[3]*Aic
Ai_3[5][2] = rho*u[3]*u[1]
Ai_3[5][3] = Aid + rho*u[3]*u[2]
Ai_3[5][4] = rho*u[3]*u[3]
Ai_3[5][5] = rho*u[3]*Aie
*********************************************************************
K_ij[15][15] = {0}

//i=1, j=1
K_ij[2][2]   = 2*mu + lambda
K_ij[3][3]   = mu
K_ij[4][4]   = mu
K_ij[5][2]   = (2*mu + lambda)*u[1]
K_ij[5][3]   = mu*u[2]
K_ij[5][4]   = mu*u[3]
K_ij[5][5]   = k


//i=1, j=2
K_ij[2][8]   = lambda
K_ij[3][7]   = mu
K_ij[5][7]   = mu*u[2]
K_ij[5][8]   = mu*u[1]


//i=1, j=3
K_ij[2][14]  = lambda
K_ij[4][12]  = mu
K_ij[5][12]  = mu*u[3]
K_ij[5][14]  = mu*u[1]


//i=2, j=1
K_ij[7][3]   = mu
K_ij[8][2]   = lambda
K_ij[10][2]  = lambda*u[2]
K_ij[10][3]  = mu*u[1]


//i=2, j=2
K_ij[7][7]   = mu
K_ij[8][8]   = 2*mu + lambda
K_ij[9][9]   = mu
K_ij[10][7]  = mu*u[1]
K_ij[10][8]  = (2*mu + lambda)*u[2]
K_ij[10][9]  = mu*u[3]
K_ij[10][10] = k


//i=2, j=3
K_ij[8][14]  = lambda
K_ij[9][13]  = mu
K_ij[10][13] = mu*u[3]
K_ij[10][14] = lambda*u[2]


//i=3, j=1
K_ij[12][4]  = mu
K_ij[14][2]  = lambda
K_ij[15][2]  = lambda*u[3]
K_ij[15][4]  = mu*u[1]


//i=3, j=2
K_ij[13][9]  = mu
K_ij[14][8]  = lambda
K_ij[15][8]  = lambda*u[3]
K_ij[15][9]  = mu*u[2]


//i=3, j=3
K_ij[12][12] = mu
K_ij[13][13] = mu
K_ij[14][14] = 2*mu + lambda
K_ij[15][12] = mu*u[1]
K_ij[15][13] = mu*u[2]
K_ij[15][14] = (2*mu + lambda)*u[3]
K_ij[15][15] = k





````









