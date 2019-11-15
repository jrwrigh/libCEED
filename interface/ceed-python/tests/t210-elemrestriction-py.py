# @file
# Test creation and view of an element restriction

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  ne = 3

  ind = np.zeros(2*ne, dtype="int32")
  for i in range(ne):
    ind[2*i+0] = i+0
    ind[2*i+1] = i+1
  r = ceed.ElemRestriction(ne, 2, ne+1, 1, ind, cmode=libceed.USE_POINTER)

  print(r)
