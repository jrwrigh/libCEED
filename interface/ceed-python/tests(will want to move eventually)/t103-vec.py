# @file
# Test setting one vector from array of another vector

import sys
from ceed import mem_host, use_pointer, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.ceed(sys.argv[1])

  n = 10

  x = ceed.vector(n)
  y = ceed.vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.setArray(mem_host, use_pointer, a)

  x_array = x.getArray(mem_host)
  y.setArray(mem_host, use_pointer, x_array)
  x.restoreArray()

  y_array = y.getArrayRead(mem_host)
  for i in range(n):
    if y_array[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array y[%d] = %f"%(i, y_array[i]))

  # LCOV_EXCL_STOP
  y.restoreArrayRead()
