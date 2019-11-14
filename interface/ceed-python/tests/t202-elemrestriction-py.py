# @file
# Test creation, use, and destruction of a blocked element restriction

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  ne = 8
  blksize = 5

  x = ceed.Vector(ne+1)
  a = np.arange(10, 10 + ne+1, dtype="float64")
  x.set_array(libceed.MEM_HOST, libceed.USE_POINTER, a)

  ind = np.zeros(2*ne, dtype="int32")
  for i in range(ne):
    ind[2*i+0] = i
    ind[2*i+1] = i+1
  r = ceed.BlockedElemRestriction(ne, 2, blksize, ne+1, 1, libceed.MEM_HOST, libceed.USE_POINTER, ind)

  y = ceed.Vector(2*blksize*2)
  y.set_value(0)

  r.apply(libceed.NOTRANSPOSE, libceed.NOTRANSPOSE, x, y, libceed.REQUEST_IMMEDIATE)

  print(y)

  x.set_value(0)
  r.apply(libceed.TRANSPOSE, libceed.NOTRANSPOSE, y, x, libceed.REQUEST_IMMEDIATE)
  print(x)
