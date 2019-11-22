# @file
# Test grad transpose with a 2D Simplex non-tensor H1 basis

import sys, math
import libceed
import numpy as np
import buildmats as bm

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  P, Q, dim = 6, 4, 2

  qref = np.empty(dim*Q, dtype="float64")
  qweight = np.empty(Q, dtype="float64")
  colsum = np.empty(P, dtype="float64")

  interp, grad = bm.buildmats(qref, qweight)

  for i in range(P):
    colsum[i] = 0
    for j in range(Q*dim):
      colsum[i] += grad[i+j*P]

  b = ceed.BasisH1(libceed.TRIANGLE, 1, P, Q, interp, grad, qref, qweight)

  in_vec = ceed.Vector(Q*dim)
  in_vec.set_value(1)
  out_vec = ceed.Vector(P)
  out_vec.set_value(0)

  b.T.apply(1, libceed.EVAL_GRAD, in_vec, out_vec)

  # Check values at quadrature points
  out_array = out_vec.get_array_read()
  for i in range(P):
    if math.fabs(colsum[i] - out_array[i]) > 1E-14:
      # LCOV_EXCL_START
      print("[%d] %f != %f"%(i, out[i], colsum[i]))
  # LCOV_EXCL_STOP
  out_vec.restore_array_read()
