# @file
# Test creation and distruction of a H1Lagrange basis

import sys
import libceed

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  b = ceed.BasisTensorH1Lagrange(1, 1, 4, 4, libceed.GAUSS_LOBATTO)
  print(b)
  del b

  b = ceed.BasisTensorH1Lagrange(1, 1, 4, 4, libceed.GAUSS)
  print(b)
  del b
