# @file
# Test creation, setting, reading, restoring, and destroying of a vector using mem_device

import sys
from libceed import MEM_HOST, MEM_DEVICE, USE_POINTER, COPY_VALUES
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)
  y = ceed.Vector(n)
  a = np.zeros(n, dtype="float64")
  x.set_array(MEM_HOST, USE_POINTER, a)

  b = x.get_array_read(MEM_DEVICE)
  y.set_array(MEM_DEVICE, COPY_VALUES, b)

  c = y.get_array_read(MEM_HOST)
  for i in range(n):
    if c[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array c[%d] = %f"%(i, c[i]))
  # LCOV_EXCL_STOP

  y.restore_array_read()
