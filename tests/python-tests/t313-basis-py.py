# @file
# Test interpolation in multiple dimensions

import sys, math
import libceed
import numpy as np

def eval(dim, x):
  result, center = 1, 0.1
  for d in range(dim):
    result *= math.tanh(x[d] - center)
    center += 0.1
  return result

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  for dim in range(1, 4):
    Q = 10
    Qdim = Q**dim
    Xdim = 2**dim
    x = np.empty(Xdim*dim, dtype="float64")
    uq = np.empty(Qdim, dtype="float64")

    for d in range(dim):
      for i in range(Xdim):
        x[d*Xdim + i] = 1 if (i % (2**(dim-d))) // (2**(dim-d-1)) else -1

    X = ceed.Vector(Xdim*dim)
    X.set_array(x, cmode=libceed.USE_POINTER)
    Xq = ceed.Vector(Qdim*dim)
    Xq.set_value(0)
    U = ceed.Vector(Qdim)
    U.set_value(0)
    Uq = ceed.Vector(Qdim)

    bxl = ceed.BasisTensorH1Lagrange(dim, dim, 2, Q, libceed.GAUSS_LOBATTO)
    bul = ceed.BasisTensorH1Lagrange(dim, 1, Q, Q, libceed.GAUSS_LOBATTO)

    bxl.apply(1, libceed.EVAL_INTERP, X, Xq)

    xq = Xq.get_array_read()
    for i in range(Qdim):
      xx = np.empty(dim, dtype="float64")
      for d in range(dim):
        xx[d] = xq[d*Qdim + i]
      uq[i] = eval(dim, xx)

    Xq.restore_array_read()
    Uq.set_array(uq, cmode=libceed.USE_POINTER)

    # This operation is the identity because the quadrature is collocated
    bul.T.apply(1, libceed.EVAL_INTERP, Uq, U)

    bxg = ceed.BasisTensorH1Lagrange(dim, dim, 2, Q, libceed.GAUSS)
    bug = ceed.BasisTensorH1Lagrange(dim, 1, Q, Q, libceed.GAUSS)

    bxg.apply(1, libceed.EVAL_INTERP, X, Xq)
    bug.apply(1, libceed.EVAL_INTERP, U, Uq)

    xq = Xq.get_array_read()
    u = Uq.get_array_read()

    for i in range(Qdim):
      xx = np.empty(dim, dtype="float64")
      for d in range(dim):
        xx[d] = xq[d*Qdim + i]
      fx = eval(dim, xx)
      if math.fabs(u[i] - fx) > 1E-4:
        # LCOV_EXCL_START
        print("[%d] %f != %f=f(%f"%(dim, u[i], fx, xx[0]))
        for d in range(dim):
          print(",%f"%xx[d])
        print(")")
        # LCOV_EXCL_START

    Xq.restore_array_read()
    Uq.restore_array_read()
