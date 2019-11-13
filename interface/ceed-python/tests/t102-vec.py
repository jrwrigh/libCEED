# @file
# Test getArrayRead state counter

import sys
from libceed import CEED_MEM_HOST
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  # Two read accesses should not generate an error
  a = x.get_array_read(CEED_MEM_HOST)
  b = x.get_array_read(CEED_MEM_HOST)

  x.restore_array_read()
  x.restore_array_read()
