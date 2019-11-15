# @file
# Test Test Symmetric Schur Decomposition

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  A = np.array([0.19996678, 0.0745459, -0.07448852, 0.0332866,
                0.0745459, 1., 0.16666509, -0.07448852,
                -0.07448852, 0.16666509, 1., 0.0745459,
                0.0332866, -0.07448852, 0.0745459, 0.19996678], dtype="float64")
  lam = np.empty(4, dtype="float64")

  libceed.Basis.symmetric_schur_decomposition(ceed, A, 4)
