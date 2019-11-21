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

from _ceed_cffi import ffi, lib
import sys
import io
from abc import ABC
from ceed_vector import Vector
from ceed_basis import BasisTensorH1, BasisTensorH1Lagrange, BasisH1
from ceed_elemrestriction import ElemRestriction, IdentityElemRestriction, BlockedElemRestriction
from ceed_qfunction import QFunction, QFunctionByName, IdentityQFunction
from ceed_operator import Operator, CompositeOperator
from ceed_constants import *

# ------------------------------------------------------------------------------
class Ceed():
  """Ceed: core components."""
  # Attributes
  _pointer = ffi.NULL

  # Constructor
  def __init__(self, resource = "/cpu/self"):
    # libCEED object
    self._pointer = ffi.new("Ceed *")

    # libCEED call
    resourceAscii = ffi.new("char[]", resource.encode("ascii"))
    lib.CeedInit(resourceAscii, self._pointer)

  # Representation
  def __repr__(self):
    return "<Ceed instance at " + hex(id(self)) + ">"

  # Get Resource
  def get_resource(self):
    """Get the full resource name for a Ceed context.

       Returns:
         resource: resource name"""

    # libCEED call
    resource = ffi.new("char **")
    lib.CeedGetResource(self._pointer[0], resource)

    return ffi.string(resource[0]).decode("UTF-8")

  # Preferred MemType
  def get_preferred_memtype(self):
    """Return Ceed preferred memory type.

       Returns:
         memtype: Ceed preferred memory type"""

    # libCEED call
    memtype = ffi.new("CeedMemType *", MEM_HOST)
    lib.CeedGetPreferredMemType(self._pointer[0], memtype)

    return memtype[0]

  # CeedVector
  def Vector(self, size):
    """Ceed Vector: storing and manipulating vectors.

       Args:
         size: length of vector

       Returns:
         vector: Ceed Vector"""

    return Vector(self, size)

  # CeedElemRestriction
  def ElemRestriction(self, nelem, elemsize, nnodes, ncomp, indices,
                      memtype=lib.CEED_MEM_HOST, cmode=lib.CEED_COPY_VALUES):
    """Ceed ElemRestriction: restriction from vectors to elements.

       Args:
         nelem: number of elements describen in the indices array
         elemsize: size (number of nodes) per element
         nnodes: the number of nodes in the local vector. The input Ceed Vector
                   to which the restriction will be applied is of size
                   nnodes * ncomp. This size may include data
                   used by other Ceed ElemRestriction objects describing
                   different types of elements.
         ncomp: number of field components per interpolation node
         indices: Numpy array of shape [nelem, elemsize]. Row i holds the
                      ordered list of the indices (into the input Ceed Vector)
                      for the unknowns corresponding to element i, where
                      0 <= i < nelem. All indices must be in the range
                      [0, nnodes - 1].
         **memtype: memory type of the indices array, default CEED_MEM_HOST
         **cmode: copy mode for the indices array, default CEED_COPY_VALUES

       Returns:
         elemrestriction: Ceed ElemRestiction"""

    return ElemRestriction(self, nelem, elemsize, nnodes, ncomp, indices,
                           memtype=memtype, cmode=cmode)

  def IdentityElemRestriction(self, nelem, elemsize, nnodes, ncomp):
    """Ceed Identity ElemRestriction: identity restriction from vectors to elements.

       Args:
         nelem: number of elements describen in the indices array
         elemsize: size (number of nodes) per element
         nnodes: the number of nodes in the local vector. The input Ceed Vector
                   to which the restriction will be applied is of size
                   nnodes * ncomp. This size may include data
                   used by other Ceed ElemRestriction objects describing
                   different types of elements.
         ncomp: number of field components per interpolation node

       Returns:
         elemrestriction: Ceed Identity ElemRestiction"""

    return IdentityElemRestriction(self, nelem, elemsize, nnodes, ncomp)

  def BlockedElemRestriction(self, nelem, elemsize, blksize, nnodes, ncomp,
                             indices, memtype=lib.CEED_MEM_HOST,
                             cmode=lib.CEED_COPY_VALUES):
    """Ceed Blocked ElemRestriction: blocked restriction from vectors to elements.

       Args:
         nelem: number of elements describen in the indices array
         elemsize: size (number of nodes) per element
         blksize: number of elements in a block
         nnodes: the number of nodes in the local vector. The input Ceed Vector
                   to which the restriction will be applied is of size
                   nnodes * ncomp. This size may include data
                   used by other Ceed ElemRestriction objects describing
                   different types of elements.
         ncomp: number of field components per interpolation node
         indices: Numpy array of shape [nelem, elemsize]. Row i holds the
                      ordered list of the indices (into the input Ceed Vector)
                      for the unknowns corresponding to element i, where
                      0 <= i < nelem. All indices must be in the range
                      [0, nnodes - 1].
         **memtype: memory type of the indices array, default CEED_MEM_HOST
         **cmode: copy mode for the indices array, default CEED_COPY_VALUES

       Returns:
         elemrestriction: Ceed Blocked ElemRestiction"""

    return BlockedElemRestriction(self, nelem, elemsize, blksize, nnodes,
                                  ncomp, indices, memtype=memtype,
                                  cmode=cmode)

  # CeedBasis
  def BasisTensorH1(self, dim, ncomp, P1d, Q1d, interp1d, grad1d,
                    qref1d, qweight1d):
    """Ceed Tensor H1 Basis: fully discrete finite element-like objects with a tensor product H^1 descretizations.

       Args:
         dim: topological dimension
         ncomp: number of field components (1 for scalar fields)
         P1d: number of nodes in one dimension
         Q1d: number of quadrature points in one dimension
         interp1d: Numpy array holding row-major Q1d × P1d matrix expressing the
                     values of nodal basis functions at quadrature points
         grad1d: Numpy array holding row-major Q1d × P1d matrix expressing the
                   derivatives of nodal basis functions at quadrature points
         qref1d: Array of length Q1d holding the locations of quadrature points
                   on the 1D reference element [-1, 1]
         qweight1d: Array of length Q1d holding the quadrature weights on the
                      reference element

       Returns:
         basis: Ceed Basis"""

    return BasisTensorH1(self, dim, ncomp, P1d, Q1d, interp1d, grad1d,
                         qref1d, qweight1d)

  def BasisTensorH1Lagrange(self, dim, ncomp, P, Q, qmode):
    """Ceed Tensor H1 Lagrange Basis: fully discrete finite element-like objects with a tensor product Lagrange basis.

       Args:
         dim: topological dimension
         ncomp: number of field components (1 for scalar fields)
         P: number of Gauss-Lobatto nodes in one dimension.  The
              polynomial degree of the resulting Q_k element is k=P-1.
         Q: number of quadrature points in one dimension
         qmode: distribution of the Q quadrature points (affects order of
                  accuracy for the quadrature)

       Returns:
         basis: Ceed Basis"""

    return BasisTensorH1Lagrange(self, dim, ncomp, P, Q, qmode)

  def BasisH1(self, topo, ncomp, nnodes, nqpts, interp, grad, qref, qweight):
    """Ceed H1 Basis: fully discrete finite element-like objects with a H^1 descretization.

       Args:
         topo: topology of the element, e.g. hypercube, simplex, etc
         ncomp: number of field components (1 for scalar fields)
         nnodes: total number of nodes
         nquts: total number of quadrature
         interp: Numpy array holding row-major nqpts x nnodes matrix expressing
                   the values of nodal basis functions at quadrature points
         grad: Numpy array holding row-major (nqpts x dim) x nnodes matrix
                 expressing the derivatives of nodal basis functions at
                 quadrature points
         qref: Array of length nqpts x dim holding the locations of quadrature
                 points on the reference element [-1, 1]
         qweight: Array of length nnodes holding the quadrature weights on the
                    reference element

       Returns:
         basis: Ceed Basis"""

    return BasisH1(self, topo, ncomp, nnodes, nqpts, interp, grad, qref, qweight)

  # CeedQFunction
  def QFunction(self, vlength, f, source):
    """Ceed QFunction: independent operations at quadrature points.

       Args:
         vlength: vector length. Caller must ensure that number of quadrature
                    points is a multiple of vlength
         f: ctypes function pointer to evaluate action at quadrature points
         source: absolute path to source of QFunction,
           "\\abs_path\\file.h:function_name

       Returns:
         qfunction: Ceed QFunction"""

    return QFunction(self, vlength, f, source)

  def QFunctionByName(self, name):
    """Ceed QFunction By Name: independent operations at quadrature points from gallery.

       Args:
         name: name of QFunction to use from gallery

       Returns:
         qfunction: Ceed QFunction By Name"""

    return QFunctionByName(self, name)

  def IdentityQFunction(self, size, inmode, outmode):
    """Ceed Idenity QFunction: identity qfunction operation.

       Args:
         size: size of the qfunction fields
         inmode: CeedEvalMode for input to Ceed QFunction
         outmode: CeedEvalMode for output to Ceed QFunction

       Returns:
         qfunction: Ceed Identity QFunction"""

    return IdentityQFunction(self, size, inmode, outmode)

  # CeedOperator
  def Operator(self, qf, dqf=None, qdfT=None):
    """Ceed Operator: composed FE-type operations on vectors.

       Args:
         qf: QFunction defining the action of the operator at quadrature points
         **dqf: QFunction defining the action of the Jacobian, default None
         **dqfT: QFunction defining the action of the transpose of the Jacobian,
                   default None

       Returns:
         operator: Ceed Operator"""

    return Operator(self, qf, dqf, qdfT)

  def CompositeOperator(self):
    """Composite Ceed Operator: composition of multiple CeedOperators.

       Returns:
         operator: Ceed Composite Operator"""

    return CompositeOperator(self)

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedDestroy(self._pointer)

# ------------------------------------------------------------------------------
