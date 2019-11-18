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
from ceed_vector import Vector
from ceed_basis import Basis, BasisTensorH1, BasisTensorH1Lagrange, BasisH1
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
    """Get the full resource name for a Ceed context."""
    # libCEED call
    resource = ffi.new("char **")
    lib.CeedGetResource(self._pointer[0], resource)

    return ffi.string(resource[0]).decode("UTF-8")

  # Preferred MemType
  def get_preferred_memtype(self):
    """Return Ceed preferred memory type."""
    # libCEED call
    memtype = ffi.new("CeedMemType *", MEM_HOST)
    lib.CeedGetPreferredMemType(self._pointer[0], memtype)

    return memtype[0]

  # CeedVector
  def Vector(self, size):
    """CeedVector: storing and manipulating vectors."""
    return Vector(self, size)

  # CeedElemRestriction
  def ElemRestriction(self, nelem, elemsize, nnodes, ncomp,indices,
                      memtype=lib.CEED_MEM_HOST, cmode=lib.CEED_COPY_VALUES):
    """CeedElemRestriction: restriction from vectors to elements."""
    return ElemRestriction(self, nelem, elemsize, nnodes, ncomp, indices,
                           memtype=memtype, cmode=cmode)

  def IdentityElemRestriction(self, nelem, elemsize, nnodes, ncomp):
    """CeedElemRestriction: identity restriction from vectors to elements."""
    return IdentityElemRestriction(self, nelem, elemsize, nnodes, ncomp)

  def BlockedElemRestriction(self, nelem, elemsize, blksize, nnodes, ncomp,
                             indices, memtype=lib.CEED_MEM_HOST,
                             cmode=lib.CEED_COPY_VALUES):
    """CeedElemRestriction: blocked restriction from vectors to elements."""
    return BlockedElemRestriction(self, nelem, elemsize, blksize, nnodes,
                                  ncomp, indices, memtype=memtype,
                                  cmode=cmode)

  # CeedBasis
  def BasisTensorH1(self, dim, ncomp, P1d, Q1d, interp1d, grad1d,
                    qref1d, qweight1d):
    """Tensor product basis class for H^1 discretizations."""
    return BasisTensorH1(self, dim, ncomp, P1d, Q1d, interp1d, grad1d,
                         qref1d, qweight1d)

  def BasisTensorH1Lagrange(self, dim, ncomp, P, Q, qmode):
    """Tensor product Lagrange basis class."""
    return BasisTensorH1Lagrange(self, dim, ncomp, P, Q, qmode)

  def BasisH1(self, topo, ncomp, nnodes, nqpts, interp, grad, qref, qweight):
    """Non tensor product basis class for H^1 discretizations."""
    return BasisH1(self, topo, ncomp, nnodes, nqpts, interp, grad, qref, qweight)

  # CeedQFunction
  def QFunction(self, vlength, f, source):
    """CeedQFunction: independent operations at quadrature points."""
    return QFunction(self, vlength, f, source)

  def QFunctionByName(self, name):
    """CeedQFunctionByName: independent operations at quadrature points from gallery."""
    return QFunctionByName(self, name)

  def IdentityQFunction(self, size, inmode, outmode):
    """CeedIdenityQFunction: identity qfunction operation."""
    return IdentityQFunction(self, size, inmode, outmode)

  # CeedOperator
  def Operator(self, qf, dqf=None, qdfT=None):
    """CeedOperator: composed FE-type operations on vectors."""
    return Operator(self, qf, dqf, qdfT)

  def CompositeOperator(self):
    """CompositeCeedOperator: composition of multiple CeedOperators."""
    return CompositeOperator(self)

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedDestroy(self._pointer)

# ------------------------------------------------------------------------------
