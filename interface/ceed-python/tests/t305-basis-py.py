# @file
# Test Simultaneous Diagonalization

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  M = np.array([0.19996678, 0.0745459, -0.07448852, 0.0332866,
                0.0745459, 1., 0.16666509, -0.07448852,
                -0.07448852, 0.16666509, 1., 0.0745459,
                0.0332866, -0.07448852, 0.0745459, 0.19996678], dtype="float64")
  K = np.array([3.03344425, -3.41501767, 0.49824435, -0.11667092,
                -3.41501767, 5.83354662, -2.9167733, 0.49824435,
                0.49824435, -2.9167733, 5.83354662, -3.41501767,
                -0.11667092, 0.49824435, -3.41501767, 3.03344425], dtype="float64")

  x, lam = libceed.Basis.simultaneous_diagonalization(ceed, K, M, 4)

  print("x: ")
  for i in range(4):
    for j in range(4):
      if x[j+4*i] <= 1E-14 and x[j+4*i] >= -1E-14:
        x[j+4*i] = 0
      print("%12.8f"%x[j+4*i])

  print("lambda: ")
  for i in range(4):
    if lam[i] <= 1E-14 and lam[i] >= -1E-14:
      lam[i] = 0
    print("%12.8f"%lam[i])
