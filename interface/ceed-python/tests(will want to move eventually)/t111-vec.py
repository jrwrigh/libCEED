# @file
# Test view

import sys
from ceed import mem_host, use_pointer, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.ceed(sys.argv[1])

  n = 10
  x = ceed.vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.setArray(mem_host, use_pointer, a)

  x.view()
