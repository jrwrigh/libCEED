# @file
# Test creation, use, and destruction of an identity element restriction

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  ne = 3

  x = ceed.Vector(2*ne)
  a = np.arange(10, 10 + 2*ne, dtype="float64")
  x.set_array(libceed.MEM_HOST, libceed.USE_POINTER, a)

  r = ceed.IdentityElemRestriction(ne, 2, 2*ne, 1)

  y = ceed.Vector(2*ne)
  y.set_value(0)

  r.apply(libceed.NOTRANSPOSE, libceed.NOTRANSPOSE, x, y, libceed.REQUEST_IMMEDIATE)

  y_array = y.get_array(libceed.MEM_HOST)
  for i in range(2*ne):
    if (10+i != y_array[i]):
      # LCOV_EXCL_START
      print("Error in restricted array y[%d] = %f"%(i, y_array[i]))
  # LCOV_EXCL_STOP
  y.restore_array()
