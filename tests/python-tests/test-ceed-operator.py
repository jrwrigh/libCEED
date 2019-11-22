# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

# @file
# Test Ceed Operator functionality

import os
import glob
import ctypes
import libceed
import numpy as np

#-------------------------------------------------------------------------------
# Utility
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Test creation creation, action, and destruction for mass matrix operator
#-------------------------------------------------------------------------------
def test_500(ceed_resource):
  ceed = libceed.Ceed(ceed_resource)

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
  file_dir = os.path.dirname(os.path.abspath(__file__))
  qfs = load_qfs_so()

  qf_setup = ceed.QFunction(1, qfs.setup_mass,
                            os.path.join(file_dir, "test-qfunctions.h:setup_mass"))
  qf_setup.add_input("weights", 1, libceed.EVAL_WEIGHT)
  qf_setup.add_input("dx", 1, libceed.EVAL_GRAD)
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

  # Apply mass matrix
  u.set_value(0)
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  for i in range(q):
    assert abs(v_array[i]) < 1E-14

  v.restore_array_read()

#-------------------------------------------------------------------------------
# Test creation creation, action, and destruction for mass matrix operator
#-------------------------------------------------------------------------------
def test_401(ceed_resource):
  ceed = libceed.Ceed(ceed_resource)

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

  qf_setup = ceed.QFunction(1, qfs.setup_mass,
                            os.path.join(file_dir, "test-qfunctions.h:setup_mass"))
  qf_setup.add_input("weights", 1, libceed.EVAL_WEIGHT)
  qf_setup.add_input("dx", 1, libceed.EVAL_GRAD)
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

  # Apply mass matrix
  u.set_value(1.)
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  total = 0.0
  for i in range(nu):
    total = total + v_array[i]
  assert abs(total - 1.0) < 1E-14

  v.restore_array_read()

#-------------------------------------------------------------------------------
