# @file
# Test sync

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

  b = np.zeros(n)
  y.set_array(libceed.MEM_HOST, libceed.USE_POINTER, b)

  c = x.get_array_read(libceed.MEM_DEVICE)
  y.set_array(libceed.MEM_DEVICE, libceed.COPY_VALUES, c)

  y.sync_array(libceed.MEM_HOST)
  for i in range(n):
    if b[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %f"%(i, b[i]))
  # LCOV_EXCL_STOP
