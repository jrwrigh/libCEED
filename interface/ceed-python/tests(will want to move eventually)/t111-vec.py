# @file
# Test view

import sys
from libceed import MEM_HOST, USE_POINTER
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(MEM_HOST, USE_POINTER, a)

  x.view()
