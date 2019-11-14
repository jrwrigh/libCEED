# @file
# Test setting one vector from array of another vector

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10

  x = ceed.Vector(n)
  y = ceed.Vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(libceed.MEM_HOST, libceed.USE_POINTER, a)

  x_array = x.get_array(libceed.MEM_HOST)
  y.set_array(libceed.MEM_HOST, libceed.USE_POINTER, x_array)
  x.restore_array()

  y_array = y.get_array_read(libceed.MEM_HOST)
  for i in range(n):
    if y_array[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array y[%d] = %f"%(i, y_array[i]))
  # LCOV_EXCL_STOP

  y.restore_array_read()
