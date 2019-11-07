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

class OperatorBase(ABC):
  # Apply CeedOperator
  def apply(self, u, v, request):
    """Apply CeedOperator to a vector."""
    # libCEED call
    libceed.CeedOperatorApply(self.op, u.vec, v.vec, request)

  # Destructor
  def __del__(self):
    # libCEED call
    libceed.CeedOperatorDestroy(self.op)

class Operator(OperatorBase):
  """CeedOperator: composed FE-type operations on vectors."""
  # Attributes
  self.ceed = fii.NULL
  self.op = ffi.NULL
  self.qf = ffi.NULL
  self.dqf = ffi.NULL
  self.dqfT = ffi.NULL
  self.fields = []

  # Constructor
  def __init__(self, ceed, qf, dqf, qdfT):
    # CeedOperator object
    self.op = ffi.new("CeedOperator *")

    # References to dependencies
    self.ceed = ceed
    self.qf = qf
    if (dqf):
      self.dqf = dqf
    if (dqfT):
      self.dqfT = dqfT

    # libCEED call
    libceed.CeedOperatorCreate(self.ceed, self.qf, self.dqf, self.dqfT, self.op)

  # References to field components
  class OperatorField:
    def __init(self, restriction, basis, vec):
      # References to dependencies
      self.rstr = restriction
      if basis:
        self.basis = basis
      if vec:
        self.vec = vec

  # Add field to CeedOperator
  def setField(self, fieldname, restriction, lmode, basis, vec):
    """Provide a field to a CeedOperator for use by its CeedQFunction."""
    # References to dependencies
    self.fields.append(OperatorField(restriction, basis, vec))

    # libCEED call
    libceed.CeedOperatorSetField(self.op, fieldname, restriction.rstr,
                                     lmode, basis.basis, vec.vec)

class CompositeOperator(OperatorBase):
  """CompositeCeedOperator: composition of multiple CeedOperators."""
  # Attributes
  self.ceed = fii.NULL
  self.op = ffi.NULL

  # Constructor
  def __init__(self, ceed):
    # CeedOperator object
    self.op = ffi.new("CeedOperator *")

    # References to dependencies
    self.ceed = ceed

    # libCEED call
    libceed.CeedCompositeOperatorCreate(self.ceed, self.op)
