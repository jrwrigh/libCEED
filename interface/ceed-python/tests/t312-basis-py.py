# @file
# Test polynomial interpolation in 1D

import sys, math
import libceed
import numpy as np

def poly_eval(x, n, p):
  y = p[n-1]
  for i in reversed(range(n-1)):
    y = y*x + p[i]

  return y

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  Q = 6
  p = [1, 2, 3, 4, 5, 6] # 1 + 2x + 3x^2 + ...
  u = np.empty(Q, dtype="float64")
  x = np.empty(2, dtype="float64")

  X = ceed.Vector(2)
  Xq = ceed.Vector(Q)
  Xq.set_value(0)
  U = ceed.Vector(Q)
  Uq = ceed.Vector(Q)
  Uq.set_value(0)
  W = ceed.Vector(Q)
  W.set_value(0)

  bxl = ceed.BasisTensorH1Lagrange(1, 1, 2, Q, libceed.GAUSS_LOBATTO)

  for i in range(len(x)):
    x[i] = (-1)**(i+1)

  X.set_array(x, cmode=libceed.USE_POINTER)
  bxl.apply(1, libceed.NOTRANSPOSE, libceed.EVAL_INTERP, X, Xq)

  xq = Xq.get_array_read()
  n = len(p)
  for i in range(Q):
    u[i] = poly_eval(xq[i], n, p)
  Xq.restore_array_read()
  U.set_array(u, cmode=libceed.USE_POINTER)

  bxg = ceed.BasisTensorH1Lagrange(1, 1, 2, Q, libceed.GAUSS)
  bug = ceed.BasisTensorH1Lagrange(1, 1, Q, Q, libceed.GAUSS)

  bxg.apply(1, libceed.NOTRANSPOSE, libceed.EVAL_INTERP, X, Xq)
  bug.apply(1, libceed.NOTRANSPOSE, libceed.EVAL_INTERP, U, Uq)
  bug.apply(1, libceed.NOTRANSPOSE, libceed.EVAL_WEIGHT, libceed.VECTOR_NONE, W)

  w = W.get_array_read()
  uq = Uq.get_array_read()
  sum = 0
  for i in range(Q):
    sum += w[i] * uq[i]
  W.restore_array_read()
  Uq.restore_array_read()

  pint = np.zeros(n+1, dtype="float64")
  for i in range(int(n)):
    pint[i+1] = p[i] / (i+1)

  error = sum - poly_eval(1, len(pint), pint) + poly_eval(-1, len(pint), pint)
  if error > 1.E-10:
    # LCOV_EXCL_START
    print("Error %e  sum %g  exact %g"%(error, sum,
          poly_eval(1, len(pint), pint) - poly_eval(-1, len(pint), pint)))
 # LCOV_EXCL_STOP
