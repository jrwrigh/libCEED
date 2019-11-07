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

from abc import ABC

class QFunctionBase(ABC):
  # Apply CeedQFunction
   def apply(self, q, u, v):
    # libCEED call
     libceed.CeedQFunctionApply(self.qf, q, u.vec, v.vec)

  # Destructor
  def __del__(self):
    # libCEED call
    libceed.CeedQFunctionDestroy(self.qf)

class QFunction(QFunctionBase):
  def __init__(self, ceed, vlength, f, source):
    # libCEED object
    self.qf = ffi.new("CeedQFunction *")

    # References to dependencies
    self.ceed = ceed

   # libCEED call
   libceed.CeedQFunctionCreateInterior(self.ceed, vlength, f, source, self.qf)

  # Add fields to CeedQFunction
  def addInput(self, fieldname, size, emode):
    # libCEED call
    libceed.CeedQFunctionAddInput(self.qf, fieldname, size, emode)

  def addOutput(self, fieldname, size, emode):
    # libCEED call
    libceed.CeedQFunctionAddOutput(self.qf, fieldname, size, emode)

class QFunctionByName(QFunctionBase):
  def __init__(self, ceed, name):
    # libCEED object
    self.qf = ffi.new("CeedQFunction *")

    # References to dependencies
    self.ceed = ceed

   # libCEED call
   libceed.CeedQFunctionCreateByName(self.ceed, name, self.qf)

class QFunctionIdentity(QFunctionBase):
  def __init__(self, ceed, size):
    # libCEED object
    self.qf = ffi.new("CeedQFunction *")

    # References to dependencies
    self.ceed = ceed

   # libCEED call
   libceed.CeedQFunctionCreateIdentity(self.ceed, size, self.qf)
