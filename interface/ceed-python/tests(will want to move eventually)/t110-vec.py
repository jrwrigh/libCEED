# @file
# Test CeedVectorGetArray state counter

import sys
from ceed import mem_host, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.ceed(sys.argv[1])

  n = 10
  x = ceed.vector(n)

  # Two write accesses should generate an error
  a = x.getArray(mem_host)
  b = x.getArray(mem_host)

  # LCOV_EXCL_START
  x.restoreArray()
  x.restoreArray()
  # LCOV_EXCL_STOP
