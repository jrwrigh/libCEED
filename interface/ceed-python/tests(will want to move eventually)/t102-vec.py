# @file
# Test getArrayRead state counter

import sys
from libceed import MEM_HOST
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  # Two read accesses should not generate an error
  a = x.GetArrayRead(MEM_HOST)
  b = x.GetArrayRead(MEM_HOST)

  x.RestoreArrayRead()
  x.RestoreArrayRead()
