# @file
# Test GetNumNodes and GetNumQuadraturePoints for basis

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  b = ceed.BasisTensorH1Lagrange( 3, 1, 4, 5, libceed.GAUSS_LOBATTO)

  P = libceed.Basis.get_num_nodes(b)
  Q = libceed.Basis.get_num_quadrature_points(b)

  if P != 64:
    # LCOV_EXCL_START
    print("%d != 64"%P)
  # LCOV_EXCL_STOP
  if Q != 125:
    # LCOV_EXCL_START
    print("%d != 125"%Q)
  # LCOV_EXCL_STOP
