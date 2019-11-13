# @file
# Test creation, setting, reading, restoring, and destroying of a vector using mem_device

import sys
from libceed import CEED_MEM_HOST, CEED_MEM_DEVICE, CEED_USE_POINTER, CEED_COPY_VALUES
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)
  y = ceed.Vector(n)
  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(CEED_MEM_HOST, CEED_USE_POINTER, a)

  b = x.get_array_read(CEED_MEM_DEVICE)
  y.set_array(CEED_MEM_DEVICE, CEED_COPY_VALUES, b)

  c = y.get_array_read(CEED_MEM_HOST)
  for i in range(n):
    if c[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array c[%d] = %f"%(i, c[i]))
  # LCOV_EXCL_STOP

  y.restore_array_read()
