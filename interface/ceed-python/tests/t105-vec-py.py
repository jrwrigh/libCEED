# @file
# Test creation, setting, reading, restoring, and destroying of a vector using mem_device

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)
  y = ceed.Vector(n)
  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(a, cmode=libceed.USE_POINTER)

  b = x.get_array_read(memtype=libceed.MEM_DEVICE)
  y.set_array(b, memtype=libceed.MEM_DEVICE, cmode=libceed.COPY_VALUES)
  x.restore_array_read(memtype=libceed.MEM_DEVICE)

  c = y.get_array_read()
  for i in range(n):
    if c[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array c[%d] = %f"%(i, c[i]))
  # LCOV_EXCL_STOP

  y.restore_array_read()
