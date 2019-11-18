# @file
# Test that length of BasisApply input/output vectors is incompatible with basis dimensions

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  Q, P, dim, ncomp = 8, 2, 3, 1
  length = Q**dim
  U = ceed.Vector(length)
  V = ceed.Vector(length + 1)

  b = ceed.BasisTensorH1Lagrange(dim, ncomp, P, Q, libceed.GAUSS)

  # Basis apply will error because dimensions don't agree
  b.apply(1, libceed.EVAL_INTERP, U, V)
