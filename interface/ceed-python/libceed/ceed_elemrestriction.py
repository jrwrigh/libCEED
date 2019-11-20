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
import numpy as np
import sys
import io
import copy
from abc import ABC
from ceed_constants import REQUEST_IMMEDIATE, REQUEST_ORDERED, MEM_HOST, COPY_VALUES, TRANSPOSE, NOTRANSPOSE
from ceed_vector import _VectorWrap

# ------------------------------------------------------------------------------
class _ElemRestrictionBase(ABC):
  """Ceed ElemRestriction: restriction from local vectors to elements."""

  # Attributes
  _ceed = ffi.NULL
  _pointer = ffi.NULL

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedElemRestrictionDestroy(self._pointer)

  # Representation
  def __repr__(self):
    return "<CeedElemRestriction instance at " + hex(id(self)) + ">"

  # String conversion for print() to stdout
  def __str__(self):
    """View an ElemRestriction via print()."""

    # libCEED call
    lib.CeedElemRestrictionView(self._pointer[0], sys.stdout)
    return ""

  # Apply CeedElemRestriction
  def apply(self, u, v, request=REQUEST_IMMEDIATE, lmode=NOTRANSPOSE,
            tmode=NOTRANSPOSE):
    """Restrict an L-vector to an E-vector or apply its transpose."""
    # libCEED call
    lib.CeedElemRestrictionApply(self._pointer[0], tmode, lmode,
                                 u._pointer[0], v._pointer[0], request)

  # Transpose an ElemRestriction
  @property
  def T(self):
    """Transpose an ElemRestriction."""

    return TransposeElemRestriction(self)

  # Transpose an ElemRestriction
  @property
  def transpose(self):
    """Transpose an ElemRestriction."""

    return TransposeElemRestriction(self)

  # Create restriction vectors
  def create_vector(self, createLvec = True, createEvec = True):
    """Create Vectors associated with an ElemRestriction."""
    # Vector pointers
    lvecPointer = ffi.new("CeedVector *") if createLvec else ffi.NULL
    evecPointer = ffi.new("CeedVector *") if createEvec else ffi.NULL

    # libCEED call
    lib.CeedElemRestrictionCreateVector(self._pointer[0], lvecPointer,
                                        evecPointer)

    # Return vectors
    lvec = _VectorWrap(self._ceed._pointer, lvecPointer) if createLvec else None
    evec = _VectorWrap(self._ceed._pointer, evecPointer) if createEvec else None

    # Return
    return [lvec, evec]

  # Get ElemRestriction multiplicity
  def get_multiplicity(self):
    """Get the multiplicity of nodes in an ElemRestriction."""
    # Create mult vector
    [mult, evec] = self.create_vector(createEvec = False)
    mult.set_value(0)

    # libCEED call
    lib.CeedElemRestrictionGetMultiplicity(self._pointer[0], mult._pointer[0])

    # Return
    return mult

# ------------------------------------------------------------------------------
class ElemRestriction(_ElemRestrictionBase):
  """Ceed ElemRestriction: restriction from local vectors to elements."""

  # Constructor
  def __init__(self, ceed, nelem, elemsize, nnodes, ncomp, indices,
               memtype=MEM_HOST, cmode=COPY_VALUES):
    # CeedVector object
    self._pointer = ffi.new("CeedElemRestriction *")

    # Reference to Ceed
    self._ceed = ceed

    # Setup the numpy array for the libCEED call
    indices_pointer = ffi.new("const CeedInt *")
    indices_pointer = ffi.cast("const CeedInt *",
                               indices.__array_interface__['data'][0])

    # libCEED call
    lib.CeedElemRestrictionCreate(self._ceed._pointer[0], nelem, elemsize,
                                  nnodes, ncomp, memtype, cmode,
                                  indices_pointer, self._pointer)

# ------------------------------------------------------------------------------
class IdentityElemRestriction(_ElemRestrictionBase):
  """Ceed IdentityElemRestriction: identity restriction from vectors to elements."""

  # Constructor
  def __init__(self, ceed, nelem, elemsize, nnodes, ncomp):
    # CeedVector object
    self._pointer = ffi.new("CeedElemRestriction *")

    # Reference to Ceed
    self._ceed = ceed

    # libCEED call
    lib.CeedElemRestrictionCreateIdentity(self._ceed._pointer[0], nelem,
                                          elemsize, nnodes, ncomp,
                                          self._pointer)

# ------------------------------------------------------------------------------
class BlockedElemRestriction(_ElemRestrictionBase):
  """Ceed Blocked ElemRestriction: blocked restriction from vectors to elements."""

  # Constructor
  def __init__(self, ceed, nelem, elemsize, blksize, nnodes, ncomp, indices,
               memtype=MEM_HOST, cmode=COPY_VALUES):
    # CeedVector object
    self._pointer = ffi.new("CeedElemRestriction *")

    # Reference to Ceed
    self._ceed = ceed

    # Setup the numpy array for the libCEED call
    indices_pointer = ffi.new("const CeedInt *")
    indices_pointer = ffi.cast("const CeedInt *",
                               indices.__array_interface__['data'][0])

    # libCEED call
    lib.CeedElemRestrictionCreateBlocked(self._ceed._pointer[0], nelem,
                                         elemsize, blksize, nnodes, ncomp,
                                         memtype, cmode, indices_pointer,
                                         self._pointer)

  # Transpose a Blocked ElemRestriction
  @property
  def T(self):
    """Transpose an ElemRestriction."""

    return TransposeBlockedElemRestriction(self)

  # Transpose a Blocked ElemRestriction
  @property
  def transpose(self):
    """Transpose an ElemRestriction."""

    return TransposeBlockedElemRestriction(self)

  # Apply CeedElemRestriction to single block
  def apply_block(self, block, u, v, lmode=NOTRANSPOSE, tmode=NOTRANSPOSE,
                  request=REQUEST_IMMEDIATE):
    """Restrict an L-vector to a block of an E-vector or apply its transpose."""
    # libCEED call
    lib.CeedElemRestrictionApplyBlock(self._pointer[0], block, tmode,
                                      lmode, u._pointer[0], v._pointer[0],
                                      request)

# ------------------------------------------------------------------------------
class TransposeElemRestriction():
  """Ceed ElemRestriction: restriction from local vectors to elements."""

  # Attributes
  _elemrestriction = None

  # Constructor
  def __init__(self, elemrestriction):

    # Reference elemrestriction
    self._elemrestriction = elemrestriction

  # Representation
  def __repr__(self):
    return "<Transpose CeedElemRestriction instance at " + hex(id(self)) + ">"


  # Apply Transpose CeedElemRestriction
  def apply(self, u, v, request=REQUEST_IMMEDIATE, lmode=NOTRANSPOSE):
    """Restrict an E-vector to an L-vector."""

    # libCEED call
    self._elemrestriction.apply(u, v, request=request, lmode=lmode,
                                tmode=TRANSPOSE)

# ------------------------------------------------------------------------------
class TransposeBlockedElemRestriction(TransposeElemRestriction):
  """Transpose Ceed Blocked ElemRestriction: blocked restriction from vectors to elements."""

  # Apply Transpose CeedElemRestriction
  def apply_block(self, block, u, v, request=REQUEST_IMMEDIATE,
                  lmode=NOTRANSPOSE):
    """Restrict a block of an E-vector to a L-vector."""

    # libCEED call
    self._elemrestriction.apply_block(block, u, v, request=request, lmode=lmode,
                                tmode=TRANSPOSE)

# ------------------------------------------------------------------------------
