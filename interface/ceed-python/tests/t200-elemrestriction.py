# @file
# Test getArray to modify array

import sys
from libceed import CEED_MEM_HOST, CEED_USE_POINTER, CEED_NOTRANSPOSE, CEED_REQUEST_IMMEDIATE
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  ne = 3

  x = ceed.Vector(ne+1)
  a = np.arange(10, 10 + ne+1, dtype="float64")
  x.set_array(CEED_MEM_HOST, CEED_USE_POINTER, a)

  ind = np.zeros(2*ne)
  for i in range(ne):
    ind[2*i+0] = i
    ind[2*i+1] = i+1
  r = ceed.ElemRestriction(ne, 2, ne+1, 1, CEED_MEM_HOST, CEED_USE_POINTER, ind)

  y = ceed.Vector(2*ne)
  y.set_value(0)

  r.apply(CEED_NOTRANSPOSE, CEED_NOTRANSPOSE, x, y, CEED_REQUEST_IMMEDIATE)

  y_array = y.get_array(MEM_HOST)
  for i in range(2*ne):
    if (10+(i+1)/2 != y_array[i]):
      # LCOV_EXCL_START
      print("Error in restricted array y[%d] = %f"%(i, y_array[i]))
  # LCOV_EXCL_STOP
  y.restore_array()
