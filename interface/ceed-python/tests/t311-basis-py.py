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
  uq = np.empty(Q, dtype="float64")
  x = np.empty(2, dtype="float64")

  X = ceed.Vector(2)
  Xq = ceed.Vector(Q)
  Xq.set_value(0)
  U = ceed.Vector(Q)
  U.set_value(0)
  Uq = ceed.Vector(Q)

  bxl = ceed.BasisTensorH1Lagrange(1, 1, 2, Q, libceed.GAUSS_LOBATTO)
  bul = ceed.BasisTensorH1Lagrange(1, 1, Q, Q, libceed.GAUSS_LOBATTO)

  for i in range(len(x)):
    x[i] = (-1)**(i+1)

  X.set_array(x, cmode=libceed.USE_POINTER)
  bxl.apply(1, libceed.EVAL_INTERP, X, Xq)

  xq = Xq.get_array_read()
  n = len(p)
  for i in range(Q):
    uq[i] = poly_eval(xq[i], n, p)
  Xq.restore_array_read()
  Uq.set_array(uq, cmode=libceed.USE_POINTER)

  # This operation is the identity because the quadrature is collocated
  bul.T.apply(1, libceed.EVAL_INTERP, Uq, U)

  bxg = ceed.BasisTensorH1Lagrange(1, 1, 2, Q, libceed.GAUSS)
  bug = ceed.BasisTensorH1Lagrange(1, 1, Q, Q, libceed.GAUSS)

  bxg.apply(1, libceed.EVAL_INTERP, X, Xq)
  bug.apply(1, libceed.EVAL_INTERP, U, Uq)

  xq = Xq.get_array_read()
  uuq = Uq.get_array_read()
  for i in range(Q):
    px = poly_eval(xq[i], n, p)
    if math.fabs(uuq[i] - px) > 1E-14:
      # LCOV_EXCL_START
      print("%f != %f=p(%f)"%(uuq[i], px, xq[i]))
    # LCOV_EXCL_STOP
  Xq.restore_array_read()
  Uq.restore_array_read()

