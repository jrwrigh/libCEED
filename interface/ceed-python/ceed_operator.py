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

from _ceed import ffi, lib
from abc import ABC

# ------------------------------------------------------------------------------
class _OperatorBase(ABC):
  # Attributes
  ceed = ffi.NULL
  pointer = ffi.NULL

  # Apply CeedOperator
  def apply(self, u, v, request):
    """Apply CeedOperator to a vector."""
    # libCEED call
    lib.CeedOperatorApply(self.pointer, u.pointer[0], v.pointer[0], request)

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedOperatorDestroy(self.pointer)

# ------------------------------------------------------------------------------
class _Operator(_OperatorBase):
  """CeedOperator: composed FE-type operations on vectors."""
  qf = ffi.NULL
  dqf = ffi.NULL
  dqfT = ffi.NULL
  fields = []

  # Constructor
  def __init__(self, ceed, qf, dqf = None, qdfT = None):
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
    lib.CeedOperatorCreate(self.ceed.pointer[0], self.qf.pointer[0],
                           self.dqf.pointer[0], self.dqfT.pointer[0],
                           self.pointer)

  # Add field to CeedOperator
  def setField(self, fieldname, restriction, lmode, basis, vector):
    """Provide a field to a CeedOperator for use by its CeedQFunction."""

    # libCEED call
    fieldnameAscii = ffi.new("char[]", fieldname.encode('ascii'))
    lib.CeedOperatorSetField(self.pointe[0], fieldnameAscii,
                             restriction.pointer[0], lmode, basis.pointer[0],
                             vector.pointer[0])

# ------------------------------------------------------------------------------
class _CompositeOperator(_OperatorBase):
  """CompositeCeedOperator: composition of multiple CeedOperators."""
  # Attributes
  subs = []

  # Constructor
  def __init__(self, ceed):
    # CeedOperator object
    self.op = ffi.new("CeedOperator *")

    # References to dependencies
    self.ceed = ceed

    # libCEED call
    lib.CeedCompositeOperatorCreate(self.ceed.pointer[0], self.pointer)

  # Add sub operators
  def addSub(self, subop):
    """Add a sub-operator to a composite CeedOperator."""
    # References to dependencies
    self.subs.append(subop)

    # libCEED call
    lib.CeedOperatorAddSup(self.pointer[0], subop.pointer[0])

# ------------------------------------------------------------------------------
