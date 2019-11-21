# @file
# Test creation and destruction of a 2D Simplex non-tensor H1 basis

import sys
import libceed
import numpy as np
import buildmats as bm

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  P, Q, dim = 6, 4, 2

  qref = np.empty(dim*Q, dtype="float64")
  qweight = np.empty(Q, dtype="float64")

  interp, grad = bm.buildmats(qref, qweight)

  b = ceed.BasisH1(libceed.TRIANGLE, 1, P, Q, interp, grad, qref, qweight)
  print(b)
