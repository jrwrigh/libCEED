# @file
# Test getArrayRead state counter

import sys
from ceed import mem_host, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  # Two read accesses should not generate an error
  a = x.GetArrayRead(mem_host)
  b = x.GetArrayRead(mem_host)

  x.RestoreArrayRead()
  x.RestoreArrayRead()
