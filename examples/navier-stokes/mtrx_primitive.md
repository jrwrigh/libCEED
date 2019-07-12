````
uu = u[1]*u[1] + u[2]*u[2] + u[3]*u[3]




A0[5][5] = {0}

A0(1, 1) = 1/(Rd*T)
A0(1, 5) = -rho/(T)

A0(2, 1) = u[1]/(Rd*T)
A0(2, 3) = rho
A0(2, 5) = -rho*u[1]/T

A0(3, 1) = u[2]/(Rd*T)
A0(3, 3) = rho
A0(3, 5) = -rho*u[2]/T

A0(4, 1) = u[3]/(Rd*T)
A0(4, 4) = rho
A0(4, 5) = -rho*u[3]/T

A0(5, 1) = (cv + uu/(2*T) + g*z/T - (cp*theta0*exp(N*N*z/g))/T)/Rd
A0(5, 2) = rho*u[1]
A0(5, 3) = rho*u[2]
A0(5, 4) = rho*u[3]
A0(5, 5) = rho*(-g*z/T - uu/(2*T) + g*gcexp(-N*N*z/g)/(theta0*N*N))
*********************************************************************

Aia = (g*z)/(Rd*T) - rho*cp*theta0
Aib = -rho*g*z/T + rho*g*g*exp(-N*N*z/g)/(theta0*N*N)
Aic = 1 + cv/Rd + g*z/Rd + uu/(2*Rd*T) - rho*cp*theta0*exp(N*N*z/g)
Aid = Pi*(1 + cv/Rd) + rho*uu/2 + rho*g*z
Aie = -g*z/T - uu/(2*T) + g*g*exp(-N*N*z/g)/(theta0*N*N)



Ai1[5][5] = {0}

Ai1(1, 1) = u[1]/(Rd*T)
Ai1(1, 2) = rho
Ai1(1, 5) = -rho*u[1]/T

Ai1(2, 1) = Aia + u[1]*u[1]/(Rd*T) + 1
Ai1(2, 2) = 2*rho*u[1]
Ai1(2, 5) = Aib - rho*u[1]*u[1]/T

Ai1(3, 1) = Aia + u[1]*u[2]/(Rd*T)
Ai1(3, 2) = rho*u[2]
Ai1(3, 3) = rho*u[1]
Ai1(2, 5) = Aib - rho*u[1]*u[2]/T

Ai1(4, 1) = Aia + u[1]*u[3]/(Rd*T)
Ai1(4, 2) = rho*u[3]
Ai1(4, 4) = rho*u[1]
Ai1(4, 5) = Aib - rho*u[1]*u[3]/T

Ai1(5, 1) = u[1]*Aic
Ai1(5, 2) = Aid + rho*u[1]*u[1]
Ai1(5, 3) = rho*u[1]*u[2]
Ai1(5, 4) = rho*u[1]*u[3]
Ai1(5, 5) = rho*u[1]*Aie




Ai2[5][5] = {0}

Ai2(1, 1) = u[2]/(Rd*T)
Ai2(1, 3) = rho
Ai2(1, 5) = -rho*u[2]/T

Ai2(2, 1) = Aia + u[2]*u[1]/(Rd*T)
Ai2(2, 2) = rho*u[2]
Ai2(2, 3) = rho*u[1]
Ai2(2, 5) = Aib - rho*u[2]*u[1]/T

Ai2(3, 1) = Aia + u[2]*u[2]/(Rd*T) +1
Ai2(3, 3) = 2*rho*u[2]
Ai2(3, 5) = Aib - rho*u[2]*u[2]/T

Ai2(4, 1) = Aia + u[2]*u[3]/(Rd*T)
Ai2(4, 3) = rho*u[3]
Ai2(4, 4) = rho*u[2]
Ai2(4, 5) = Aib - rho*u[2]*u[3]/T

Ai2(5, 1) = u[2]*Aic
Ai2(5, 2) = rho*u[2]*u[1]
Ai2(5, 3) = Aid + rho*u[2]*u[2]
Ai2(5, 4) = rho*u[2]*u[3]
Ai2(5, 5) = rho*u[2]*Aie



Ai3[5][5] = {0}

Ai3(1, 1) = u[3]/(Rd*T)
Ai3(1, 4) = rho
Ai3(1, 5) = -rho*u[3]/T

Ai3(2, 1) = Aia + u[3]*u[1]/(Rd*T)
Ai3(2, 2) = rho*u[3]
Ai3(2, 4) = rho*u[1]
Ai3(2, 5) = Aib - rho*u[3]*u[1]/T

Ai3(3, 1) = Aia + u[3]*u[2]/(Rd*T)
Ai3(3, 3) = rho*u[3]
Ai3(3, 4) = rho*u[2]
Ai3(3, 5) = Aib - rho*u[3]*u[2]/T

Ai3(4, 1) = Aia + u[3]*u[3]/(Rd*T) +1
Ai3(4, 4) = 2*rho*u[3]
Ai3(4, 5) = Aib - rho*u[3]*u[3]/T

Ai3(5, 1) = u[3]*Aic
Ai3(5, 2) = rho*u[3]*u[1]
Ai3(5, 3) = Aid + rho*u[3]*u[2]
Ai3(5, 4) = rho*u[3]*u[3]
Ai3(5, 5) = rho*u[3]*Aie

````









