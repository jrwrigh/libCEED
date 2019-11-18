# @file
# Test viewing of qfunction

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

  print(qf_setup)
  print(qf_mass)
