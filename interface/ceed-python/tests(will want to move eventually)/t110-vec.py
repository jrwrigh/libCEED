# @file
# Test CeedVectorGetArray state counter

import sys
from libceed import MEM_HOST
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  # Two write accesses should generate an error
  a = x.get_array(MEM_HOST)
  b = x.get_array(MEM_HOST)

  # LCOV_EXCL_START
  x.restore_array()
  x.restore_array()
  # LCOV_EXCL_STOP
