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
from abc import ABC

# ------------------------------------------------------------------------------
class _ElemRestrictionBase(ABC):

  # Attributes
  _ceed = ffi.NULL
  _pointer = ffi.NULL

  # Apply CeedElemRestriction
  def apply(self, tmode, lmode, u, v, request):
    """Restrict an L-vector to an E-vector or apply its transpose."""
    # libCEED call
    lib.CeedElemRestrictionApply(self._pointer, tmode, lmode, u._pointer[0],
                                 v._pointer[0], request)

  # Create restriction vectors
  def create_vector(self, createLvec = True, createEvec = True):
    """Create CeedVectors associated with a ElemRestriction."""
    # Vector pointers
    lvecPointer = ffi.new("CeedVector *") if createLvec else ffi.NULL
    evecPointer = ffi.new("CeedVector *") if createEvec else ffi.NULL

    # libCEED call
    lib.CeedElemRestrictionCreateVector(self._pointer[0], lvecPointer,
                                        evecPointer)

    # Return vectors
    lvec = _VectorClone(self._ceed[0], lvecPointer) if createLvec else None
    evec = _VectorClone(self._ceed[0], evecPointer) if createEvec else None

    # Return
    return [lvec, evec]

  # Get ElemRestriction multiplicity
  def get_multiplicity(self):
    """Get the multiplicity of nodes in a ElemRestriction."""
    # Create mult vector
    [mult, evec] = self.createVector(createEvec = False)

    # libCEED call
    lib.CeedElemRestrictionApply(self._pointer[0], mult._pointer[0])

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedElemRestrictionDestroy(self._pointer)

# ------------------------------------------------------------------------------
class ElemRestriction(_ElemRestrictionBase):
  """Ceed ElemRestriction: restriction from vectors to elements."""

  # Constructor
  def __init__(self, ceed, nelem, elemsize, nnodes, ncomp, mtype, cmode,
               indices):
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
                                  nnodes, ncomp, mtype, cmode, indices_pointer,
                                  self._pointer)

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
  """Ceed BlockedElemRestriction: blocked restriction from vectors to elements."""

  # Constructor
  def __init__(self, ceed, nelem, elemsize, blksize, nnodes, ncomp, mtype,
               cmode, indices):
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
                                         mtype, cmode, indices_pointer,
                                         self._pointer)

  # Apply CeedElemRestriction to single block
  def apply_block(self, block, tmode, lmode, u, v, request):
    """Restrict an L-vector to a block of an E-vector or apply its transpose."""
    # libCEED call
    lib.CeedElemRestrictionApplyBlock(self._pointer, block, tmode, lmode,
                                      u._pointer[0], v._pointer[0], request)


# ------------------------------------------------------------------------------
