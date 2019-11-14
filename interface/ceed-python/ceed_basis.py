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
import numpy as np
from abc import ABC

# ------------------------------------------------------------------------------
class _BasisBase(ABC):
  """Ceed Basis: class for creating and manipulating finite element bases."""

  # Attributes
  _ceed = ffi.NULL
  _pointer = ffi.NULL

  # Representation
  def __repr__(self):
    return "<CeedBasis instance at " + hex(id(self)) + ">"

  # String conversion for print() to stdout
  def __str__(self):
    """View a Basis via print()."""

    # libCEED call
    lib.CeedBasisView(self._pointer[0], sys.stdout)
    return ""

  # Apply Basis
  def apply(self, nelem, tmode, emode, u, v):
    """Apply basis evaluation from nodes to quadrature points or viceversa"""
    # libCEED call
    CeedBasisApply(self._pointer[0], nelem, tmode, emode,
                   u._pointer[0], v._pointer[0])

  # Gauss quadrature
  @staticmethod
  def gauss_quadrature(Q):
    """Construct a Gauss-Legendre quadrature"""

    # Setup numpy arrays
    qref1d = np.empty(Q, dtype="float64")
    qweight1d = np.empy(Q, dtype="float64")

    # Setup pointers
    qref1d_pointer = ffi.new("CeedScalar *")
    qref1d_pointer = ffi.cast("CeedScalar *", qref1d.__array_interface__['data'][0])
    qweight1d_pointer = ffi.new("CeedScalar *")
    qweight1d_pointer = ffi.cast("CeedScalar *", qweight1d.__array_interface__['data'][0])

    # libCEED call
    lib.CeedGaussQuadrature(Q, qref1d_pointer, qweight1d_pointer)

    return qref1d, qweight1d

  # Lobatto quadrature
  @staticmethod
  def lobatto_quadrature(Q):
    """Construct a Gauss-Legendre-Lobatto quadrature"""

    # Setup arguments
    qref1d = np.empty(Q, dtype="float64")
    qref1d_pointer = ffi.new("CeedScalar *")
    qref1d_pointer = ffi.cast("CeedScalar *", qref1d.__array_interface__['data'][0])
    qweight1d = np.empy(Q, dtype="float64")
    qweight1d_pointer = ffi.new("CeedScalar *")
    qweight1d_pointer = ffi.cast("CeedScalar *", qweight1d.__array_interface__['data'][0])

    # libCEED call
    lib.CeedLobattoQuadrature(Q, qref1d_pointer, qweight1d_pointer)

    return qref1d, qweight1d

  # Scalar view
  @staticmethod
  def scalar_view(name, m, n, array, format = ffi.NUL, file = sys.stdout):
    """View an array stored in a CeedBasis"""

    # Check if format is a string before encoding it
    if type(format) == "str":
      fstr = format.encode("ascii", "strict")
    else:
      fstr = format

    # Setup arguments
    a_pointer = ffi.new("CeedScalar *")
    a_pointer = ffi.cast("CeedScalar *", array.__array_interface__['data'][0])
    namestr = name.encode("ascii", "strict")

    # libCEED call
    CeedScalarView(namestr, fstr, m, n, a_pointer, file)

  # Householder reflection
  @staticmethod
  def householder_reflect(A, v, b, m, n, row, col):
    """Compute Householder reflection"""

    # Setup arguments
    A_pointer = ffi.new("CeedScalar *")
    A_pointer = ffi.cast("CeedScalar *", A.__array_interface__['data'][0])
    v_pointer = ffi.new("CeedScalar *")
    v_pointer = ffi.cast("CeedScalar *", v.__array_interface__['data'][0])

    # libCEED call
    lib.CeedHouseholderReflect(A_pointer, v_pointer, b, m, n, row, col)

    return A

  # Apply Householder Q matrix
  @staticmethod
  def householder_apply_q(A, Q, tau, tmode, m, n, k, row, col):
    """Apply Householder Q matrix"""

    # Setup arguments
    A_pointer = ffi.new("CeedScalar *")
    A_pointer = ffi.cast("CeedScalar *", A.__array_interface__['data'][0])
    Q_pointer = ffi.new("CeedScalar *")
    Q_pointer = ffi.cast("CeedScalar *", Q.__array_interface__['data'][0])
    tau_pointer = ffi.new("CeedScalar *")
    tau_pointer = ffi.cast("CeedScalar *", tau.__array_interface__['data'][0])

    # libCEED call
    lib.CeedHouseholderApplyQ(A_pointer, Q_pointer, tau_pointer, tmode, m, n,
                              k, row, col)

    return A

  # Compute Givens rotation
  @staticmethod
  def givens_rotation(A, c, s, tmode, i, k, m, n):
    """Compute Givens rotation"""

    # Setup arguments
    A_pointer = ffi.new("CeedScalar *")
    A_pointer = ffi.cast("CeedScalar *", A.__array_interface__['data'][0])

    # libCEED call
    lib.CeedGivensRotation(A_pointer, c, s, tmode, i, k, m, n)

    return A

  # QR factorization
  @staticmethod
  def qr_factorization(ceed, mat, tau, m, n):
    """Return QR Factorization of matrix"""

    # Setup arguments
    mat_pointer = ffi.new("CeedScalar *")
    mat_pointer = ffi.cast("CeedScalar *", mat.__array_interface__['data'][0])
    tau_pointer = ffi.new("CeedScalar *")
    tau_pointer = ffi.cast("CeedScalar *", tau.__array_interface__['data'][0])

    # libCEED call
    lib.CeedQRFactorization(ceed._pointer[0], mat_pointer, tau_pointer, m, n)

    return mat, tau

  # Symmetric Schur decomposition
  @staticmethod
  def symmetric_schur_decomposition(ceed, mat, n):
    """Return symmetric Schur decomposition of a symmetric matrix mat"""

    # Setup arguments
    mat_pointer = ffi.new("CeedScalar *")
    mat_pointer = ffi.cast("CeedScalar *", mat.__array_interface__['data'][0])
    l = np.empty(n, dtype="float64")
    l_pointer = ffi.new("CeedScalar *")
    l_pointer = ffi.cast("CeedScalar *", l.__array_interface__['data'][0])

    # libCEED call
    lib.CeedSymmetricSchurDecomposition(ceed._pointer[0], mat_pointer, l_pointer, n)

    return l


