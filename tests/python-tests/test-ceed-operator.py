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
import buildmats as bm

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
def test_501(ceed_resource):
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
# Test creation creation, action, and destruction for mass matrix operator
#-------------------------------------------------------------------------------
def test_502(ceed_resource):
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
  u = ceed.Vector(2*nu)
  v = ceed.Vector(2*nu)

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
  ru = ceed.ElemRestriction(nelem, p, nu, 2, indu, cmode=libceed.USE_POINTER)
  rui = ceed.IdentityElemRestriction(nelem, q, q*nelem, 1)

  # Bases
  bx = ceed.BasisTensorH1Lagrange(1, 1, 2, q, libceed.GAUSS)
  bu = ceed.BasisTensorH1Lagrange(1, 2, p, q, libceed.GAUSS)

  # QFunctions
  file_dir = os.path.dirname(os.path.abspath(__file__))
  qfs = load_qfs_so()

  qf_setup = ceed.QFunction(1, qfs.setup_mass,
                            os.path.join(file_dir, "test-qfunctions.h:setup_mass"))
  qf_setup.add_input("weights", 1, libceed.EVAL_WEIGHT)
  qf_setup.add_input("dx", 1, libceed.EVAL_GRAD)
  qf_setup.add_output("rho", 1, libceed.EVAL_NONE)

  qf_mass = ceed.QFunction(1, qfs.apply_mass_two,
                           os.path.join(file_dir, "test-qfunctions.h:apply_mass_two"))
  qf_mass.add_input("rho", 1, libceed.EVAL_NONE)
  qf_mass.add_input("u", 2, libceed.EVAL_INTERP)
  qf_mass.add_output("v", 2, libceed.EVAL_INTERP)

  # Operators
  op_setup = ceed.Operator(qf_setup)
  op_setup.set_field("weights", rxi, bx, libceed.VECTOR_NONE)
  op_setup.set_field("dx", rx, bx, libceed.VECTOR_ACTIVE)
  op_setup.set_field("rho", rui, libceed.BASIS_COLLOCATED,
                     libceed.VECTOR_ACTIVE)

  op_mass = ceed.Operator(qf_mass)
  op_mass.set_field("rho", rui, libceed.BASIS_COLLOCATED, qdata)
  op_mass.set_field("u", ru, bu, libceed.VECTOR_ACTIVE, lmode=libceed.TRANSPOSE)
  op_mass.set_field("v", ru, bu, libceed.VECTOR_ACTIVE, lmode=libceed.TRANSPOSE)

  # Setup
  op_setup.apply(x, qdata)

  # Apply mass matrix
  u_array = u.get_array()
  for i in range(nu):
    u_array[2*i] = 1.
    u_array[2*i+1] = 2.
  u.restore_array()
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  total_1 = 0.0
  total_2 = 0.0
  for i in range(nu):
    total_1 = total_1 + v_array[2*i]
    total_2 = total_2 + v_array[2*i+1]
  assert abs(total_1 - 1.0) < 1E-13
  assert abs(total_2 - 2.0) < 1E-13

  v.restore_array_read()

#-------------------------------------------------------------------------------
# Test creation, action, and destruction for mass matrix operator with passive
#   inputs and outputs
#-------------------------------------------------------------------------------
def test_503(ceed_resource):
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
  op_mass.set_field("u", ru, bu, u, lmode=libceed.TRANSPOSE)
  op_mass.set_field("v", ru, bu, v, lmode=libceed.TRANSPOSE)

  # Setup
  op_setup.apply(x, qdata)

  # Apply mass matrix
  u.set_value(1)
  op_mass.apply(libceed.VECTOR_NONE, libceed.VECTOR_NONE)

  # Check
  v_array = v.get_array_read()
  total = 0.0
  for i in range(nu):
    total = total + v_array[i]
  assert abs(total - 1.0) < 1E-13

  v.restore_array_read()

#-------------------------------------------------------------------------------
# Test viewing of mass matrix operator
#-------------------------------------------------------------------------------
def test_504(ceed_resource, capsys):
  ceed = libceed.Ceed(ceed_resource)

  nelem = 15
  p = 5
  q = 8
  nx = nelem + 1
  nu = nelem*(p-1) + 1

  # Vectors
  qdata = ceed.Vector(nelem*q)

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

  # View
  print(op_setup)
  print(op_mass)

  stdout, stderr = capsys.readouterr()
  with open(os.path.abspath("./output/test_504.out")) as output_file:
    true_output = output_file.read()

  assert stdout == true_output

#-------------------------------------------------------------------------------
# Test creation creation, action, and destruction for mass matrix operator
#-------------------------------------------------------------------------------
def test_510(ceed_resource):
  ceed = libceed.Ceed(ceed_resource)

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

  # Apply mass matrix
  u.set_value(0.)
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  for i in range(ndofs):
    assert abs(v_array[i]) < 1E-14

  v.restore_array_read()

#-------------------------------------------------------------------------------
# Test creation creation, action, and destruction for mass matrix operator
#-------------------------------------------------------------------------------
def test_511(ceed_resource):
  ceed = libceed.Ceed(ceed_resource)

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

  # Apply mass matrix
  u.set_value(1.)
  op_mass.apply(u, v)

  # Check
  v_array = v.get_array_read()
  total = 0.0
  for i in range(ndofs):
    total = total + v_array[i]
  assert abs(total - 1.0) < 1E-10

  v.restore_array_read()

#-------------------------------------------------------------------------------
