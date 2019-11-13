# @file
# Test view

import sys
from ceed import mem_host, use_pointer, ceed
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  n = 10
  x = ceed.Vector(n)

  a = np.arange(10, 10 + n, dtype="float64")
  x.SetArray(mem_host, use_pointer, a)

  x.View()
