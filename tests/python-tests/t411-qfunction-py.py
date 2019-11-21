# @file
# Test creation, evaluation, and destruction of identity qfunction

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  qf = ceed.IdentityQFunction(1, libceed.EVAL_INTERP, libceed.EVAL_INTERP)

  q = 8

  u_array = np.zeros(q, dtype="float64")
  for i in range(q):
    u_array[i] = i*i

  u = ceed.Vector(q)
  u.set_array(u_array, cmode=libceed.USE_POINTER)
  v = ceed.Vector(q)
  v.set_value(0)

  inputs = [ u ]
  outputs = [ v ]
  qf.apply(q, inputs, outputs)

  v_array = v.get_array_read()
  for i in range(q):
    if v_array[i] != i*i:
      # LCOV_EXCL_START
      print("[%d] v %f != v_true %f"%(i, v_array[i], v_true[i]))
  # LCOV_EXCL_STOP
  v.restore_array_read()
