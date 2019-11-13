# @file
# Test getArray to modify array

import sys
from ceed import mem_host, use_pointer, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10

  x = ceed.Vector(n)
  a = np.zeros(n, dtype="float64")
  x.set_array(MEM_HOST, USE_POINTER, a)

  b = x.get_array(MEM_HOST)
  b[3] = -3.14;
  x.restore_array()

  if a[3] != -3.14:
    # LCOV_EXCL_START
    print("Error writing array a[3] = %f"%b[3])
  # LCOV_EXCL_STOP
