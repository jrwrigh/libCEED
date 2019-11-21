# @file
# Test square Gauss Lobatto interp1d is identity

import sys, math
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  dim, P1d, Q1d = 2, 4, 4
  length = int(Q1d**dim + 0.4)

  U = ceed.Vector(length)
  V = ceed.Vector(length)

  u = np.ones(length, dtype="float64")
  U.set_array(u, cmode=libceed.USE_POINTER)

  b = ceed.BasisTensorH1Lagrange(dim, 1, P1d, Q1d, libceed.GAUSS_LOBATTO)
  b.apply(1, libceed.EVAL_INTERP, U, V)

  v = V.get_array_read()
  for i in range(length):
    if math.fabs(v[i] - 1.) > 1E-15:
      # LCOV_EXCL_START
      print("v[%d] = %f != 1."%(i, v[i]))
  # LCOV_EXCL_STOP
  V.restore_array_read()



