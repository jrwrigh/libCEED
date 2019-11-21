# @file
# Test creation, evaluation, and destruction of identity qfunction with size>1

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  size = 3
  qf = ceed.IdentityQFunction(size, libceed.EVAL_INTERP, libceed.EVAL_INTERP)

  q = 8

  u_array = np.zeros(q*size, dtype="float64")
  for i in range(q*size):
    u_array[i] = i*i

  u = ceed.Vector(q*size)
  u.set_array(u_array, cmode=libceed.USE_POINTER)
  v = ceed.Vector(q*size)
  v.set_value(0)

  inputs = [ u ]
  outputs = [ v ]
  qf.apply(q, inputs, outputs)

  v_array = v.get_array_read()
  for i in range(q*size):
    if v_array[i] != i*i:
      # LCOV_EXCL_START
      print("[%d] v %f != v_true %f"%(i, v_array[i], v_true[i]))
  # LCOV_EXCL_STOP
  v.restore_array_read()
