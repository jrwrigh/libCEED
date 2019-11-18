# @file
# Test creation, evaluation, and destruction for qfunction

import sys
import os
import libceed
import ctypes
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  qfs = ctypes.cdll.LoadLibrary("./qfs.so")
  file_dir = os.path.abspath(__file__)

  qf_setup = ceed.QFunction(1, qfs.setup,
                            os.path.join(file_dir, "t400-qfunction.h:setup"))
  qf_setup.add_input("w", 1, libceed.EVAL_WEIGHT)
  qf_setup.add_output("qdata", 1, libceed.EVAL_INTERP)

  qf_mass = ceed.QFunction(1, qfs.mass,
                           os.path.join(file_dir, "t400-qfunction.h:mass"))
  qf_mass.add_input("qdata", 1, libceed.EVAL_NONE)
  qf_mass.add_input("u", 1, libceed.EVAL_INTERP)
  qf_mass.add_output("v", 1, libceed.EVAL_INTERP)

  q = 8

  w_array = np.zeros(q, dtype="float64")
  u_array = np.zeros(q, dtype="float64")
  v_true  = np.zeros(q, dtype="float64")
  for i in range(q):
    x = 2.*i/(q-1) - 1
    w_array[i] = 1 - x*x
    u_array[i] = 2 + 3*x + 5*x*x
    v_true[i]  = w_array[i] * u_array[i]

  w = ceed.Vector(q)
  w.set_array(w_array, cmode=libceed.USE_POINTER)
  u = ceed.Vector(q)
  u.set_array(u_array, cmode=libceed.USE_POINTER)
  v = ceed.Vector(q)
  v.set_value(0)
  qdata = ceed.Vector(q)
  qdata.set_value(0)

  inputs = [ w ]
  outputs = [ qdata ]
  qf_setup.apply(q, inputs, outputs)

  inputs = [ qdata, u ]
  outputs = [ v ]
  qf_mass.apply(q, inputs, outputs)

  v_array = v.get_array_read()
  for i in range(q):
    if v_array[i] != v_true[i]:
      # LCOV_EXCL_START
      print("[%d] v %f != v_true %f"%(i, v_array[i], v_true[i]))
  # LCOV_EXCL_STOP
  v.restore_array_read()
