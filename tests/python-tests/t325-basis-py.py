# @file
# Test grad transpose with a 2D Simplex non-tensor H1 basis

import sys, math
import libceed
import numpy as np
import buildmats as bm

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  P, Q, dim, ncomp = 6, 4, 2, 3

  qref = np.empty(dim*Q, dtype="float64")
  qweight = np.empty(Q, dtype="float64")
  colsum = np.empty(P, dtype="float64")

  interp, grad = bm.buildmats(qref, qweight)

  for i in range(P):
    colsum[i] = 0
    for j in range(Q*dim):
      colsum[i] += grad[i+j*P]

  b = ceed.BasisH1(libceed.TRIANGLE, ncomp, P, Q, interp, grad, qref, qweight)

  in_vec = ceed.Vector(Q*dim*ncomp)
  in_array = in_vec.get_array()
  for d in range(dim):
    for n in range(ncomp):
      for q in range(Q):
        in_array[q+(n+d*ncomp)*Q] = n*1.0
  in_vec.restore_array()
  out_vec = ceed.Vector(P*ncomp)
  out_vec.set_value(0)

  b.T.apply(1, libceed.EVAL_GRAD, in_vec, out_vec)

  # Check values at quadrature points
  out_array = out_vec.get_array_read()
  for p in range(P):
    for n in range(ncomp):
      if math.fabs(n*colsum[p] - out_array[p+n*P]) > 1E-14:
        # LCOV_EXCL_START
        print("[%d] %f != %f"%(p, out_array[p+n*P], n*colsum[p]))
  # LCOV_EXCL_STOP
  out_vec.restore_array_read()
