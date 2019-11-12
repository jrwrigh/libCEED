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

# ------------------------------------------------------------------------------
class _Vector:
  """CeedVector: storing and manipulating vectors."""

  # Attributes
  _ceed = ffi.NULL
  _pointer = ffi.NULL

  # Constructor
  def __init__(self, ceed, size):
    # CeedVector object
    self._pointer = ffi.new("CeedVector *")

    # Reference to Ceed
    self._ceed = ceed

    # libCEED call
    lib.CeedVectorCreate(self._ceed._pointer[0], size, self._pointer)

  # Set Vector's data array
  def setArray(self, mtype, cmode, array):
    """Set the array used by a CeedVector, freeing any previously allocated
       array if applicable."""
    # Setup the numpy array for the libCEED call
    array_pointer = ffi.new("CeedScalar *")
    array_pointer = ffi.cast("CeedScalar *", array.__array_interface__['data'][0])

    # libCEED call
    lib.CeedVectorSetArray(self._pointer[0], mtype, cmode, array_pointer)

  # Get Vector's data array
  def getArray(self, mtype, array):
    """Get read/write access to a CeedVector via the specified memory type."""
    # Setup the pointer's pointer
    array_pointer = ffi.new("CeedScalar **")

    # libCEED call
    lib.CeedVectorGetArray(self._pointer[0], mtype, array_pointer)
    array.__array_interface__['data'] = array_pointer

  # Get Vector's data array in read-only mode
  def getArrayRead(self, mtype):
    """Get read-only access to a CeedVector via the specified memory type."""
    # Setup the pointer's pointer
    array_pointer = ffi.new("CeedScalar **")

    # libCEED call
    lib.CeedVectorGetArrayRead(self._pointer[0], mtype, array_pointer)
    length_pointer = ffi.new("CeedInt *")
    lib.CeedVectorGetLength(self._pointer[0], length_pointer)
    ret = np.ctypeslib.as_array(array_pointer[0], length_pointer[0])
    return np.ctypeslib.as_array(array_pointer[0], length_pointer[0])
#    array.__array_interface__['data'] = array_pointer

  # Restore the Vector's data array
  def restoreArray(self, array):
    """Restore an array obtained using getArray()."""
    # Setup the pointer's pointer
    array_pointer = ffi.new("CeedScalar **")

    # libCEED call
    lib.CeedVectorRestoreArray(self._pointer[0], array_pointer)

  # Restore an array obtained using getArrayRead
  def restoreArrayRead(self, array):
    """Restore an array obtained using getArrayRead()."""
    # Setup the pointer's pointer
    array_pointer = ffi.new("CeedScalar **")

    # libCEED call
    lib.CeedVectorRestoreArrayRead(self._pointer[0], array_pointer)

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedVectorDestroy(self._pointer)

# ------------------------------------------------------------------------------
class _VectorClone:
  """Copy a CeedVector """

  # Constructor
  def __init__(self, ceed, pointer):
    # CeedVector object
    self._pointer = pointer

    # Reference to Ceed
    self._ceed = ceed

# ------------------------------------------------------------------------------
