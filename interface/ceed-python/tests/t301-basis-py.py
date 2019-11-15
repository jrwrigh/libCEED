# @file
# Test QR Factorization

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  qr = np.array([1, -1, 4, 1, 4, -2, 1, 4, 2, 1, -1, 0], dtype="float64")
  tau = np.empty(3, dtype="float64")

  qr, tau = libceed.Basis.qr_factorization(ceed, qr, tau, 4, 3)

  for i in range(len(qr)):
    if qr[i] <= 1E-14  and qr[i] >= -1E-14:
      qr[i] = 0
    print("%12.8f"%qr[i])

  for i in range(len(tau)):
    if tau[i] <= 1E-14  and tau[i] >= -1E-14:
      tau[i] = 0
    print("%12.8f"%tau[i])
