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
from ceed_vector import Vector
from ceed_elemrestriction import ElemRestriction, IdentityElemRestriction, BlockedElemRestriction
from ceed_qfunction import QFunction, QFunctionByName, IdentityQFunction
from ceed_operator import Operator, CompositeOperator

# ------------------------------------------------------------------------------
# Ceed Enums
# ------------------------------------------------------------------------------
# CeedMemType
CEED_MEM_HOST       = lib.CEED_MEM_HOST
CEED_MEM_DEVICE     = lib.CEED_MEM_DEVICE
ceed_memtypes  = {CEED_MEM_HOST:     "ceed_mem_host",
                  CEED_MEM_DEVICE:   "ceed_mem_device"}

# CeedCopyMode
CEED_COPY_VALUES    = lib.CEED_COPY_VALUES
CEED_USE_POINTER    = lib.CEED_USE_POINTER
CEED_OWN_POINTER    = lib.CEED_OWN_POINTER
ceed_copymodes = {CEED_COPY_VALUES:  "ceed_copy_values",
                  CEED_USE_POINTER:  "ceed_use_pointer",
                  CEED_OWN_POINTER:  "ceed_own_pointer"}

# CeedEvalMode
CEED_EVAL_NONE      = lib.CEED_EVAL_NONE
CEED_EVAL_INTERP    = lib.CEED_EVAL_INTERP
CEED_EVAL_GRAD      = lib.CEED_EVAL_GRAD
CEED_EVAL_DIV       = lib.CEED_EVAL_DIV
CEED_EVAL_CURL      = lib.CEED_EVAL_CURL
CEED_EVAL_WEIGHT    = lib.CEED_EVAL_WEIGHT
ceed_evalmodes = {CEED_EVAL_NONE:    "ceed_eval_none",
                  CEED_EVAL_INTERP:  "ceed_eval_interp",
                  CEED_EVAL_GRAD:    "ceed_eval_grad",
                  CEED_EVAL_DIV:     "ceed_eval_div",
                  CEED_EVAL_CURL:    "ceed_eval_curl",
                  CEED_EVAL_WEIGHT:  "ceed_eval_weight"}

# CeedElemTopology
CEED_LINE           = lib.CEED_LINE
CEED_TRIANGLE       = lib.CEED_TRIANGLE
CEED_QUAD           = lib.CEED_QUAD
CEED_TET            = lib.CEED_TET
CEED_PYRAMID        = lib.CEED_PYRAMID
CEED_PRISM          = lib.CEED_PRISM
CEED_HEX            = lib.CEED_HEX
ceed_elemtopologies = {CEED_LINE:     "ceed_line",
                       CEED_TRIANGLE: "ceed_triangle",
                       CEED_QUAD:     "ceed_quad",
                       CEED_TET:      "ceed_tet",
                       CEED_PYRAMID:  "ceed_pyramid",
                       CEED_PRISM:    "ceed_prism",
                       CEED_HEX:      "ceed_hex"}

# CeedTransposeMode
CEED_TRANSPOSE      = lib.CEED_TRANSPOSE
CEED_NOTRANSPOSE    = lib.CEED_NOTRANSPOSE
ceed_transposemodes = {CEED_TRANSPOSE:   "ceed_transpose",
                       CEED_NOTRANSPOSE: "ceed_notranspose"}

# ------------------------------------------------------------------------------
# Ceed Constants
# ------------------------------------------------------------------------------
# Requests
CEED_REQUEST_IMMEDIATE = lib.CEED_REQUEST_IMMEDIATE
CEED_REQUEST_ORDERED = lib.CEED_REQUEST_ORDERED

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
    memtype = ffi.new("CeedMemType *", CEED_MEM_HOST)
    lib.CeedGetPreferredMemType(self._pointer[0], memtype)

    return memtype[0]

  # CeedVector
  def Vector(self, size):
    """CeedVector: storing and manipulating vectors."""
    return Vector(self, size)

  # CeedElemRestriction
  def ElemRestriction(self, nelem, elemsize, nnodes, ncomp, mtype, cmode,
                      indices):
    """CeedElemRestriction: restriction from vectors to elements."""
    return ElemRestriction(self, nelem, elemsize, nnodes, ncomp, mtype,
                            cmode, indices)

  def IdentityElemRestriction(self, nelem, elemsize, nnodes, ncomp, mtype,
                              cmode):
    """CeedElemRestriction: identity restriction from vectors to elements."""
    return IdentityElemRestriction(self, nelem, elemsize, nnodes, ncomp, mtype,
                                    cmode)

  def BlockedElemRestriction(self, nelem, elemsize, blksize, nnodes, ncomp,
                             mtype, cmode, indices):
    """CeedElemRestriction: blocked restriction from vectors to elements."""
    return BlockedElemRestriction(self, nelem, elemsize, blksize, nnodes,
                                   ncomp, mtype, cmode, indices)

  # CeedBasis

  # CeedQFunction
  def QFunction(self, vlength, f, source):
    """CeedQFunction: independent operations at quadrature points."""
    return QFunction(self, vlength, f, source)

  def QFunctionByName(self, name):
    """CeedQFunctionByName: independent operations at quadrature points from gallery."""
    return QFunctionByName(self, name)

  def IdentityQFunction(self, size):
    """CeedIdenityQFunction: identity qfunction operation."""
    return IdentityQFunction(self, size)

  # CeedOperator
  def Operator(self, qf, dqf, qdfT):
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
