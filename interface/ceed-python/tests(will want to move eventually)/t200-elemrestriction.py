# @file
# Test getArray to modify array

import sys
from libceed import MEM_HOST, USE_POINTER, NOTRANSPOSE, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  ne = 3

  x = ceed.Vector(ne+1)
  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(MEM_HOST, USE_POINTER, a)

  ind = np.zeros(2*ne)
  for i in range(n):
    ind[2*i+0] = i
    ind[2*i+1] = i+1
  r = ceed.ElemRestrictionCreate(ne, 2, ne+1, 1, MEM_HOST, USE_POINTER, ind)

  y = ceed.Vector(2*ne)
  y.set_value(0)

  r.apply(NOTRANSPOSE, NOTRASPOSE, x, y, REQUEST_IMMEDIATE)

  y_array = y.get_array(MEM_HOST)
  for i in range(2*ne):
    if (10+(i+1)/2 != y_array[i]):
      # LCOV_EXCL_START
      print("Error in restricted array y[%d] = %f"%(i, y_array[i]))
  # LCOV_EXCL_STOP
  y.restore_array()
