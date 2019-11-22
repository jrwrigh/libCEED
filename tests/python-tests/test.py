# @file
# Test grad transpose with a 2D Simplex non-tensor H1 basis

import sys, math
import libceed
import os
import glob
import ctypes
import libceed
import numpy as np
import buildmats as bm

def load_qfs_so():
  # Filename
  file_dir = os.path.dirname(os.path.abspath(__file__))

  # Rename, if needed
  qfs_so = glob.glob("libceed_qfunctions.*.so")
  if len(qfs_so) > 0:
    os.rename(qfs_so[0], file_dir + "/qfs.so")

  # Load library
  qfs = ctypes.cdll.LoadLibrary('./qfs.so')

  return qfs

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  nelem = 12
  dim = 2
  p = 6
  q = 4
  nx = 3
  ny = 2
  ndofs = (nx*2+1)*(ny*2+1)
  nqpts = nelem*q

  # Vectors
  x = ceed.Vector(dim*ndofs)
  x_array = np.zeros(dim*ndofs)
  for i in range(ndofs):
    x_array[i] = (1. / (nx*2)) * (i % (nx*2+1))
    x_array[i+ndofs] = (1. / (ny*2)) * (i / (nx*2+1))
  x.set_array(x_array, cmode=libceed.USE_POINTER)

  qdata = ceed.Vector(nqpts)
  u = ceed.Vector(ndofs)
  v = ceed.Vector(ndofs)

  # Restrictions
  indx = np.zeros(nelem*p, dtype="int32")
  for i in range(nelem//2):
    col = i % nx;
    row = i // nx;
    offset = col*2 + row*(nx*2+1)*2

    indx[i*2*p+ 0] =  2 + offset
    indx[i*2*p+ 1] =  9 + offset
    indx[i*2*p+ 2] = 16 + offset
    indx[i*2*p+ 3] =  1 + offset
    indx[i*2*p+ 4] =  8 + offset
    indx[i*2*p+ 5] =  0 + offset

    indx[i*2*p+ 6] = 14 + offset
    indx[i*2*p+ 7] =  7 + offset
    indx[i*2*p+ 8] =  0 + offset
    indx[i*2*p+ 9] = 15 + offset
    indx[i*2*p+10] =  8 + offset
    indx[i*2*p+11] = 16 + offset

  rx = ceed.ElemRestriction(nelem, p, ndofs, dim, indx,
                            cmode=libceed.USE_POINTER)
  rxi = ceed.IdentityElemRestriction(nelem, p, nelem*p, dim)

  ru = ceed.ElemRestriction(nelem, p, ndofs, 1, indx, cmode=libceed.USE_POINTER)
  rui = ceed.IdentityElemRestriction(nelem, q, nqpts, 1)

  # Bases
  qref = np.empty(dim*q, dtype="float64")
  qweight = np.empty(q, dtype="float64")

  interp, grad = bm.buildmats(qref, qweight)
  bx = ceed.BasisH1(libceed.TRIANGLE, dim, p, q, interp, grad, qref, qweight)

  interp, grad = bm.buildmats(qref, qweight)
  bu = ceed.BasisH1(libceed.TRIANGLE, 1, p, q, interp, grad, qref, qweight)

  # QFunctions
  file_dir = os.path.dirname(os.path.abspath(__file__))
  qfs = load_qfs_so()

  qf_setup = ceed.QFunction(1, qfs.setup_mass_2d,
                            os.path.join(file_dir, "test-qfunctions.h:setup_mass_2d"))
  qf_setup.add_input("weights", 1, libceed.EVAL_WEIGHT)
  qf_setup.add_input("dx", dim*dim, libceed.EVAL_GRAD)
  qf_setup.add_output("rho", 1, libceed.EVAL_NONE)

  qf_mass = ceed.QFunction(1, qfs.apply_mass,
                           os.path.join(file_dir, "test-qfunctions.h:apply_mass"))
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
  print(qdata)
  print(x)
  # Apply mass matrix
  u.set_value(1.)
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  total = 0.0
  for i in range(ndofs):
    total = total + v_array[i]
  print( abs(total - 1.0))

  v.restore_array_read()
