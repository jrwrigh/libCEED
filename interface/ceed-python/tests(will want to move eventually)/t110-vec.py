# @file
# Test CeedVectorGetArray state counter

import sys
from ceed import mem_host, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  # Two write accesses should generate an error
  a = x.GetArray(mem_host)
  b = x.GetArray(mem_host)

  # LCOV_EXCL_START
  x.RestoreArray()
  x.RestoreArray()
  # LCOV_EXCL_STOP
