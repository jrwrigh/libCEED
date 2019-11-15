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

  lam = libceed.Basis.symmetric_schur_decomposition(ceed, A, 4)

  print("Q: ")
  for i in range(4):
    for j in range(4):
      if A[j+4*i] <= 1E-14 and A[j+4*i] >= -1E-14:
         A[j+4*i] = 0
      print("%12.8f"%A[j+4*i])

  print("lmbda: ")
  for i in range(4):
    if lam[i] <= 1E-14 and lam[i] >= -1E-14:
      lam[i] = 0
    print("%12.8f"%lam[i])
