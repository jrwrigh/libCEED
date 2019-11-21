# @file
# Test creation, setting, reading, restoring, and destroying of a vector

import sys
import libceed
import numpy as np


if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(a, cmode=libceed.USE_POINTER)

  b = x.get_array_read()
  for i in range(n):
    if b[i] != 10 + i:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %4.2f" %(i, b[i]))
  # LCOV_EXCL_STOP

  x.restore_array_read()
