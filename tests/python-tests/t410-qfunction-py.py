# @file
# Test creation, evaluation, and destruction for qfunction by name

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  qf_setup = ceed.QFunctionByName("Mass1DBuild")
  qf_mass = ceed.QFunctionByName("MassApply")

  q = 8

  j_array = np.zeros(q, dtype="float64")
  w_array = np.zeros(q, dtype="float64")
  u_array = np.zeros(q, dtype="float64")
  v_true  = np.zeros(q, dtype="float64")
  for i in range(q):
    x = 2.*i/(q-1) - 1
    j_array[i] = 1
    w_array[i] = 1 - x*x
    u_array[i] = 2 + 3*x + 5*x*x
    v_true[i]  = w_array[i] * u_array[i]

  j = ceed.Vector(q)
  j.set_array(j_array, cmode=libceed.USE_POINTER)
  w = ceed.Vector(q)
  w.set_array(w_array, cmode=libceed.USE_POINTER)
  u = ceed.Vector(q)
  u.set_array(u_array, cmode=libceed.USE_POINTER)
  v = ceed.Vector(q)
  v.set_value(0)
  qdata = ceed.Vector(q)
  qdata.set_value(0)

  inputs = [ j, w ]
  outputs = [ qdata ]
  qf_setup.apply(q, inputs, outputs)

  inputs = [ w, u ]
  outputs = [ v ]
  qf_mass.apply(q, inputs, outputs)

  v_array = v.get_array_read()
  for i in range(q):
    if v_array[i] != v_true[i]:
      # LCOV_EXCL_START
      print("[%d] v %f != v_true %f"%(i, v_array[i], v_true[i]))
  # LCOV_EXCL_STOP
  v.restore_array_read()
