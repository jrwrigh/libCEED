# @file
# Test creation, setting, reading, restoring, and destroying of a vector

import sys
from ceed import mem_host, use_pointer
import libceed
import numpy as np


if __name__ == "__main__":
  ceed = libceed.ceed(sys.argv[1])

  n = 10
  x = ceed.vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.setArray(mem_host, use_pointer, a)

  b = x.getArrayRead(mem_host)
  for i in range(n):
    if b[i] != 10 + i:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %4.2f" %(i, b[i]))
      # LCOV_EXCL_STOP

  x.restoreArrayRead()
