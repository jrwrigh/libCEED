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
  _ceed = ffi.NULL
  _pointer = ffi.NULL

  # Apply CeedOperator
  def apply(self, u, v, request):
    """Apply CeedOperator to a vector."""
    # libCEED call
    lib.CeedOperatorApply(self._pointer[0], u._pointer[0], v._pointer[0],
                          request)

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedOperatorDestroy(self._pointer)

# ------------------------------------------------------------------------------
class _Operator(_OperatorBase):
  """CeedOperator: composed FE-type operations on vectors."""

  # Constructor
  def __init__(self, ceed, qf, dqf = None, qdfT = None):
    # CeedOperator object
    self._pointer = ffi.new("CeedOperator *")

    # Reference to Ceed
    self._ceed = ceed

    # libCEED call
    lib.CeedOperatorCreate(self._ceed._pointer[0], qf._pointer[0],
                           dqf._pointer[0] if dqf else ffi.NULL,
                           dqfT._pointer[0] if dqfT else ffi.NULL,
                           self._pointer)

  # Add field to CeedOperator
  def set_field(self, fieldname, restriction, lmode, basis, vector):
    """Provide a field to a CeedOperator for use by its CeedQFunction."""

    # libCEED call
    fieldnameAscii = ffi.new("char[]", fieldname.encode('ascii'))
    lib.CeedOperatorSetField(self._pointer[0], fieldnameAscii,
                             restriction._pointer[0], lmode, basis._pointer[0],
                             vector._pointer[0])

# ------------------------------------------------------------------------------
class _CompositeOperator(_OperatorBase):
  """CompositeCeedOperator: composition of multiple CeedOperators."""

  # Constructor
  def __init__(self, ceed):
    # CeedOperator object
    self.pointer = ffi.new("CeedOperator *")

    # Reference to Ceed
    self._ceed = ceed

    # libCEED call
    lib.CeedCompositeOperatorCreate(self._ceed._pointer[0], self._pointer)

  # Add sub operators
  def add_sub(self, subop):
    """Add a sub-operator to a composite CeedOperator."""

    # libCEED call
    lib.CeedOperatorAddSup(self._pointer[0], subop._pointer[0])

# ------------------------------------------------------------------------------
