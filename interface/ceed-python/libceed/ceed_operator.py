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
import sys
import io
from abc import ABC
from ceed_constants import REQUEST_IMMEDIATE, REQUEST_ORDERED, NOTRANSPOSE

# ------------------------------------------------------------------------------
class _OperatorBase(ABC):
  """Ceed Operator: composed FE-type operations on vectors."""

  # Attributes
  _ceed = ffi.NULL
  _pointer = ffi.NULL

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedOperatorDestroy(self._pointer)

  # Representation
  def __repr__(self):
    return "<CeedOperator instance at " + hex(id(self)) + ">"

  # String conversion for print() to stdout
  def __str__(self):
    """View a Operator via print()."""

    # libCEED call
    lib.CeedOperatorView(self._pointer[0], sys.stdout)
    return ""

  # Apply CeedOperator
  def apply(self, u, v, request=REQUEST_IMMEDIATE):
    """Apply Operator to a vector."""
    # libCEED call
    lib.CeedOperatorApply(self._pointer[0], u._pointer[0], v._pointer[0],
                          request)

# ------------------------------------------------------------------------------
class Operator(_OperatorBase):
  """Ceed Operator: composed FE-type operations on vectors."""

  # Constructor
  def __init__(self, ceed, qf, dqf = None, dqfT = None):
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
  def set_field(self, fieldname, restriction, basis, vector, lmode=NOTRANSPOSE):
    """Provide a field to a Operator for use by its QFunction."""

    # libCEED call
    fieldnameAscii = ffi.new("char[]", fieldname.encode('ascii'))
    lib.CeedOperatorSetField(self._pointer[0], fieldnameAscii,
                             restriction._pointer[0], lmode, basis._pointer[0],
                             vector._pointer[0])

# ------------------------------------------------------------------------------
class CompositeOperator(_OperatorBase):
  """Ceed Composite Operator: composition of multiple Operators."""

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
