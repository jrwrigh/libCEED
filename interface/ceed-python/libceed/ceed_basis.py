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
  """Base class for manipulating finite element bases."""

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
    """Apply basis evaluation from nodes to quadrature points or viceversa."""

    # libCEED call
    lib.CeedBasisApply(self._pointer[0], nelem, tmode, emode,
                       u._pointer[0], v._pointer[0])

  # Get number of nodes
  def get_num_nodes(self):
    """Get total number of nodes (in dim dimensions) of a Basis."""

    # Setup argument
    p_pointer = ffi.new("CeedInt *")

    # libCEED call
    lib.CeedBasisGetNumNodes(self._pointer[0], p_pointer)

    return p_pointer[0]

  # Get number of quadrature points
  def get_num_quadrature_points(self):
    """Get total number of quadrature points (in dim dimensions) of a Basis"""

    # Setup argument
    q_pointer = ffi.new("CeedInt *")

    # libCEED call
    lib.CeedBasisGetNumQuadraturePoints(self._pointer[0], q_pointer)

    return q_pointer[0]

# ------------------------------------------------------------------------------

class Basis(_BasisBase):
  """Class for Basis static methods"""

  # Gauss quadrature
  @staticmethod
  def gauss_quadrature(Q):
    """Construct a Gauss-Legendre quadrature."""

    # Setup arguments
    qref1d = np.empty(Q, dtype="float64")
    qweight1d = np.empy(Q, dtype="float64")

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
    """Construct a Gauss-Legendre-Lobatto quadrature."""

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
  def scalar_view(name, m, n, array, format = ffi.NULL, file = sys.stdout):
    """View an array stored in a CeedBasis."""

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

  # QR factorization
  @staticmethod
  def qr_factorization(ceed, mat, tau, m, n):
    """Return QR Factorization of matrix."""

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
    """Return symmetric Schur decomposition of a symmetric matrix mat."""

    # Setup arguments
    mat_pointer = ffi.new("CeedScalar *")
    mat_pointer = ffi.cast("CeedScalar *", mat.__array_interface__['data'][0])

    l = np.empty(n, dtype="float64")
    l_pointer = ffi.new("CeedScalar *")
    l_pointer = ffi.cast("CeedScalar *", l.__array_interface__['data'][0])

    # libCEED call
    lib.CeedSymmetricSchurDecomposition(ceed._pointer[0], mat_pointer, l_pointer, n)

    return l

  # Simultaneous Diagonalization
  @staticmethod
  def simultaneous_diagonalization(ceed, matA, matB, n):
    """Return Simultaneous Diagonalization of two matrices."""

    # Setup arguments
    matA_pointer = ffi.new("CeedScalar *")
    matA_pointer = ffi.cast("CeedScalar *", matA.__array_interface__['data'][0])

    matB_pointer = ffi.new("CeedScalar *")
    matB_pointer = ffi.cast("CeedScalar *", matB.__array_interface__['data'][0])

    l = np.empty(n, dtype="float64")
    l_pointer = ffi.new("CeedScalar *")
    l_pointer = ffi.cast("CeedScalar *", l.__array_interface__['data'][0])

    x = np.empty(n, dtype="float64")
    x_pointer = ffi.new("CeedScalar *")
    x_pointer = ffi.cast("CeedScalar *", x.__array_interface__['data'][0])

    # libCEED call
    lib.CeedSimultaneousDiagonalization(ceed, matA_pointer, matB_pointer,
                                        x_pointer, l_pointer, n)

    return x, l

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedBasisDestroy(self._pointer)

# ------------------------------------------------------------------------------
class BasisTensorH1(Basis):
  """Tensor product basis class for H^1 discretizations."""

  # Constructor
  def __init__(self, ceed, dim, ncomp, P1d, Q1d, interp1d, grad1d,
               qref1d, qweight1d):

    # Setup arguments
    self._pointer = ffi.new("CeedBasis *")

    self._ceed = ceed

    interp1d_pointer = ffi.new("CeedScalar *")
    interp1d_pointer = ffi.cast("CeedScalar *", interp1d.__array_interface__['data'][0])

    grad1d_pointer = ffi.new("CeedScalar *")
    grad1d_pointer = ffi.cast("CeedScalar *", grad1d.__array_interface__['data'][0])

    qref1d_pointer = ffi.new("CeedScalar *")
    qref1d_pointer = ffi.cast("CeedScalar *", qref1d.__array_interface__['data'][0])

    qweight1d_pointer = ffi.new("CeedScalar *")
    qweight1d_pointer = ffi.cast("CeedScalar *", qweight1d.__array_interface__['data'][0])

    # libCEED call
    lib.CeedBasisCreateTensorH1(self._ceed._pointer[0], dim, ncomp, P1d, Q1d,
                                interp1d_pointer, grad1d_pointer, qref1d_pointer,
                                qweight1d_pointer, self._pointer)

# ------------------------------------------------------------------------------
class BasisTensorH1Lagrange(Basis):
  """Tensor product Lagrange basis class."""

  # Constructor
  def __init__(self, ceed, dim, ncomp, P, Q, qmode):

    # Setup arguments
    self._pointer = ffi.new("CeedBasis *")

    self._ceed = ceed

    # libCEED call
    lib.CeedBasisCreateTensorH1Lagrange(self._ceed._pointer[0], dim, ncomp, P,
                                        Q, qmode, self._pointer)

# ------------------------------------------------------------------------------
class BasisH1(Basis):
  """Non tensor product basis class for H^1 discretizations."""

  # Constructor
  def __init__(self, topo, ncomp, nnodes, nqpts, interp, grad, qref, qweight):

    # Setup arguments
    self._pointer = ffi.new("CeedBasis *")

    self._ceed = ceed

    interp_pointer = ffi.new("CeedScalar *")
    interp_pointer = ffi.cast("CeedScalar *", interp.__array_interface__['data'][0])

    grad_pointer = ffi.new("CeedScalar *")
    grad_pointer = ffi.cast("CeedScalar *", grad.__array_interface__['data'][0])

    qref_pointer = ffi.new("CeedScalar *")
    qref_pointer = ffi.cast("CeedScalar *", qref.__array_interface__['data'][0])

    qweight_pointer = ffi.new("CeedScalar *")
    qweight_pointer = ffi.cast("CeedScalar *", qweight.__array_interface__['data'][0])

    # libCEED call
    lib.CeedBasisCreateH1(self._ceed._pointer[0], topo, ncomp, nnodes, nqpts,
                          interp_pointer, grad_pointer, qref_pointer,
                          qweight_pointer, self._pointer)

# ------------------------------------------------------------------------------
