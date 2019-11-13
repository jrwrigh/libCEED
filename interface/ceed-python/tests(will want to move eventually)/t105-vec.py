# @file
# Test creation, setting, reading, restoring, and destroying of a vector using mem_device

import sys
from ceed import mem_host, mem_device, use_pointer, copy_values, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.ceed(sys.argv[1])

  n = 10
  x = ceed.vector(n)
  y = ceed.vector(n)
  a = np.zeros(n, dtype="float64")
  x.setArray(mem_host, use_pointer, a)

  b = x.getArrayRead(mem_device)
  y.setArray(mem_device, copy_values, b)

  c = y.getArrayRead(mem_host)
  for i in range(n):
    if c[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array c[%d] = %f"%(i, c[i]))
  # LCOV_EXCL_STOP
  y.restoreArrayRead()
