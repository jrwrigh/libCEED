# @file
# Test grad in multiple dimensions

import sys, math
import libceed
import numpy as np

def eval(dim, x):
  result = math.tanh(x[0] + 0.1)
  if dim > 1:
    result += math.atan(x[1] + 0.2)
  if dim > 2:
    result += math.exp(-(x[2] + 0.3)*(x[2] + 0.3))
  return result

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  for dim in range (1, 4):
    P, Q = 8, 10
    Pdim = P**dim
    Qdim = Q**dim
    Xdim = 2**dim
    sum1 = sum2 = 0
    x = np.empty(Xdim*dim, dtype="float64")
    u = np.empty(Pdim, dtype="float64")

    for d in range(dim):
      for i in range(Xdim):
        x[d*Xdim + i] = 1 if (i % (2**(dim-d))) // (2**(dim-d-1)) else -1

    X = ceed.Vector(Xdim*dim)
    X.set_array(x, cmode=libceed.USE_POINTER)
    Xq = ceed.Vector(Pdim*dim)
    Xq.set_value(0)
    U = ceed.Vector(Pdim)
    Uq = ceed.Vector(Qdim*dim)
    Uq.set_value(0)
    Ones = ceed.Vector(Qdim*dim)
    Ones.set_value(1)
    Gtposeones = ceed.Vector(Pdim)
    Gtposeones.set_value(0)

    # Get function values at quadrature points
    bxl = ceed.BasisTensorH1Lagrange(dim, dim, 2, P, libceed.GAUSS_LOBATTO)
    bxl.apply(1, libceed.NOTRANSPOSE, libceed.EVAL_INTERP, X, Xq)

    xq = Xq.get_array_read()
    for i in range(Pdim):
      xx = np.empty(dim, dtype="float64")
      for d in range(dim):
        xx[d] = xq[d*Pdim + i]
      u[i] = eval(dim, xx)

    Xq.restore_array_read()
    U.set_array(u, cmode=libceed.USE_POINTER)

    # Calculate G u at quadrature points, G' * 1 at dofs
    bug = ceed.BasisTensorH1Lagrange(dim, 1, P, Q, libceed.GAUSS)
    bug.apply(1, libceed.NOTRANSPOSE, libceed.EVAL_GRAD, U, Uq)
    bug.apply(1, libceed.TRANSPOSE, libceed.EVAL_GRAD, Ones, Gtposeones)

    # Check if 1' * G * u = u' * (G' * 1)
    gtposeones = Gtposeones.get_array_read()
    uq = Uq.get_array_read()

    for i in range(Pdim):
      sum1 += gtposeones[i]*u[i]
    for i in range(dim*Qdim):
      sum2 += uq[i]
    Gtposeones.restore_array_read()
    Uq.restore_array_read()

    if math.fabs(sum1 - sum2) > 1E-10:
      # LCOV_EXCL_START
      printf("[%d] %f != %f"%(dim, sum1, sum2))
    # LCOV_EXCL_STOP
