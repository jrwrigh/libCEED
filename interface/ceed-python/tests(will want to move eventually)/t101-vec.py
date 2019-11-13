# @file
# Test setValue

import sys
from libceed import MEM_HOST, USE_POINTER
import libceed
import numpy as np

def check_values(ceed, x, value):
  b = x.get_array_read(mem_host)
  for i in range(len(b)):
    if b[i] != value:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %f"%(i,b[i]))
      # LCOV_EXCL_STOP
  x.restore_array_read()

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])
  n = 10
  x = ceed.Vector(n)
  value = 1
  a = np.arange(10, 10 + n, dtype="float64")
  x.set_array(MEM_HOST, USE_POINTER, a)

  b = x.get_array_read(MEM_HOST)
  for i in range(len(b)):
    if b[i] != 10+i:
      # LCOV_EXCL_START
      print("Error reading array b[%d] = %f"%(i,b[i]))
  # LCOV_EXCL_STOP

  x.restore_array_read()

  x.set_value(3.0)
  check_values(ceed, x, 3.0)
  del x

  x = ceed.Vector(n)
  # Set value before setting or getting the array
  x.set_value(5.0)
  check_values(ceed, x, 5.0)
