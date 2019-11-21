# @file
# Test CeedVectorGetArray state counter

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  # Two write accesses should generate an error
  a = x.get_array()
  b = x.get_array()

  # LCOV_EXCL_START
  x.restore_array()
  x.restore_array()
  # LCOV_EXCL_STOP
