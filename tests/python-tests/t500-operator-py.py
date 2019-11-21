# @file
# Test creation creation, action, and destruction for mass matrix operator

import sys
import os
import libceed
import ctypes
import numpy as np

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  nelem = 15
  p = 5
  q = 8
  nx = nelem + 1
  nu = nelem*(p-1) + 1

  # Vectors
  x = ceed.Vector(nx)
  x_array = np.zeros(nx)
  for i in range(nx):
    x_array[i] = i / (nx - 1.0)
  x.set_array(x_array, cmode=libceed.USE_POINTER)

  qdata = ceed.Vector(nelem*q)
  u = ceed.Vector(nu)
  v = ceed.Vector(nu)

  # Restrictions
  indx = np.zeros(nx*2, dtype="int32")
  for i in range(nx):
    indx[2*i+0] = i
    indx[2*i+1] = i+1
  rx = ceed.ElemRestriction(nelem, 2, nx, 1, indx, cmode=libceed.USE_POINTER)
  rxi = ceed.IdentityElemRestriction(nelem, 2, nelem*2, 1)

  indu = np.zeros(nelem*p, dtype="int32")
  for i in range(nelem):
    for j in range(p):
      indu[p*i+j] = i*(p-1) + j
  ru = ceed.ElemRestriction(nelem, p, nu, 1, indu, cmode=libceed.USE_POINTER)
  rui = ceed.IdentityElemRestriction(nelem, q, q*nelem, 1)

  # Bases
  bx = ceed.BasisTensorH1Lagrange(1, 1, 2, q, libceed.GAUSS)
  bu = ceed.BasisTensorH1Lagrange(1, 1, p, q, libceed.GAUSS)

  # QFunctions
  qfs = ctypes.cdll.LoadLibrary("./qfs.so")
  file_dir = os.path.abspath(__file__)

  qf_setup = ceed.QFunction(1, qfs.setup,
                            os.path.join(file_dir, "t500-operator.h:setup"))
  qf_setup.add_input("weights", 1, libceed.EVAL_WEIGHT)
  qf_setup.add_input("dx", 1, libceed.EVAL_GRAD)
  qf_setup.add_output("rho", 1, libceed.EVAL_NONE)

  qf_mass = ceed.QFunction(1, qfs.mass,
                           os.path.join(file_dir, "t500-operator.h:mass"))
  qf_mass.add_input("rho", 1, libceed.EVAL_NONE)
  qf_mass.add_input("u", 1, libceed.EVAL_INTERP)
  qf_mass.add_output("v", 1, libceed.EVAL_INTERP)

  # Operators
  op_setup = ceed.Operator(qf_setup)
  op_setup.set_field("weights", rxi, bx, libceed.VECTOR_NONE)
  op_setup.set_field("dx", rx, bx, libceed.VECTOR_ACTIVE)
  op_setup.set_field("rho", rui, libceed.BASIS_COLLOCATED,
                     libceed.VECTOR_ACTIVE)

  op_mass = ceed.Operator(qf_mass)
  op_mass.set_field("rho", rui, libceed.BASIS_COLLOCATED, qdata)
  op_mass.set_field("u", ru, bu, libceed.VECTOR_ACTIVE)
  op_mass.set_field("v", ru, bu, libceed.VECTOR_ACTIVE)

  # Setup
  op_setup.apply(x, qdata)

  # Apply mass matrix
  u.set_value(0)
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  for i in range(q):
    if abs(v_array[i]) > 1E-14:
      # LCOV_EXCL_START
      print("[%d] v %f != 0.0"%(i, v_array[i]))
  # LCOV_EXCL_STOP
  v.restore_array_read()
