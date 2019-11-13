# @file
# Test getArrayRead state counter

import sys
from ceed import mem_host, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.ceed(sys.argv[1])

  n = 10
  x = ceed.vector(n)

  # Two read accesses should not generate an error
  a = x.getArrayRead(mem_host)
  b = x.getArrayRead(mem_host)

  x.restoreArrayRead()
  x.restoreArrayRead()
