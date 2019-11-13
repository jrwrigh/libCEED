# @file
# Test setValue

import sys
from ceed import mem_host, use_pointer, ceed
import libceed
import numpy as np

def checkValues(ceed, x, value):
  b = x.GetArrayRead(mem_host)
  for i in range(len(b)):
    if b[i] != value:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %f"%(i,b[i]))
      # LCOV_EXCL_STOP
  x.RestoreArrayRead()

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])
  n = 10
  x = ceed.Vector(n)
  value = 1
  a = np.arange(10, 10 + n, dtype="float64")
  x.SetArray(mem_host, use_pointer, a)

  b = x.GetArrayRead(mem_host)
  for i in range(len(b)):
    if b[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %f"%(i,b[i]))
      # LCOV_EXCL_STOP

  x.RestoreArrayRead()

  x.SetValue(3.0)
  checkValues(ceed, x, 3.0)
  del x

  x = ceed.Vector(n)
  # Set value before setting or getting the array
  x.SetValue(5.0)
  checkValues(ceed, x, 5.0)
