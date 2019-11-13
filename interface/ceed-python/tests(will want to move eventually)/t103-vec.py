# @file
# Test setting one vector from array of another vector

import sys
from ceed import mem_host, use_pointer, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10

  x = ceed.Vector(n)
  y = ceed.Vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.SetArray(mem_host, use_pointer, a)

  x_array = x.GetArray(mem_host)
  y.SetArray(mem_host, use_pointer, x_array)
  x.RestoreArray()

  y_array = y.GetArrayRead(mem_host)
  for i in range(n):
    if y_array[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array y[%d] = %f"%(i, y_array[i]))
  # LCOV_EXCL_STOP

  y.RestoreArrayRead()
