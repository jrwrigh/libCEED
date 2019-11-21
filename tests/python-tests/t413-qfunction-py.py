# @file
# Test viewing of qfunction by name

import sys
import libceed
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  qf_setup = ceed.QFunctionByName("Mass1DBuild")
  qf_mass = ceed.QFunctionByName("MassApply")

  print(qf_setup)
  print(qf_mass)
