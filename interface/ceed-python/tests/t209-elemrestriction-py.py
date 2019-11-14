# @file
# Test creation, use, and destruction of an element restriction

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  ne = 3

  ind = np.zeros(4*ne, dtype="int32")
  for i in range(ne):
    ind[4*i+0] = i*3+0
    ind[4*i+1] = i*3+1
    ind[4*i+2] = i*3+2
    ind[4*i+3] = i*3+3
  r = ceed.ElemRestriction(ne, 4, 3*ne+1, 1, libceed.MEM_HOST, libceed.USE_POINTER, ind)

  mult = r.get_multiplicity()

  mult_array = mult.get_array(libceed.MEM_HOST)
  for i in range(3*ne+1):
    val = 1 + (1 if (i > 0 and i < 3*ne and i%3 == 0) else 0)
    if (val != mult_array[i]):
      # LCOV_EXCL_START
      print("Error in multiplicity array mult[%d] = %f"%(i, mult_array[i]))
  # LCOV_EXCL_STOP
  mult.restore_array()
