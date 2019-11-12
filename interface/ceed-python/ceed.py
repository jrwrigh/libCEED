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
from ceed_vector import _Vector
from ceed_qfunction import _QFunction,
                           _QFunctionByName,
                           _QFunctionIdentity
from ceed_operator import _Operator,
                          _CompositeOperator

# ------------------------------------------------------------------------------
# Ceed Enums
# ------------------------------------------------------------------------------
# CeedMemType
mem_host       = lib.CEED_MEM_HOST
mem_device     = lib.CEED_MEM_DEVICE
ceed_memtypes  = {mem_host:      "ceed_mem_host",
                  mem_device:    "ceed_mem_device"}

# CeedCopyMode
copy_values    = lib.CEED_COPY_VALUES
use_pointer    = lib.CEED_USE_POINTER
own_pointer    = lib.CEED_OWN_POINTER
ceed_memtypes  = {copy_values:    "ceed_copy_values",
                  use_pointer:    "ceed_use_pointer",
                  own_pointer:    "ceed_own_pointer"}

# CeedTransposeMode
transpose      = lib.CEED_TRANSPOSE
notranspose    = lib.CEED_NOTRANSPOSE
ceed_memtypes  = {transpose:      "ceed_transpose",
                  notranspose:    "ceed_notranspose"}

# CeedEvalMode
eval_none      = lib.CEED_EVAL_NONE
eval_interp    = lib.CEED_EVAL_INTERP
eval_grad      = lib.CEED_EVAL_GRAD
eval_div       = lib.CEED_EVAL_DIV
eval_curl      = lib.CEED_EVAL_CURL
eval_weight    = lib.CEED_EVAL_WEIGHT
ceed_memtypes  = {eval_none:      "ceed_eval_none",
                  eval_interp:    "ceed_eval_interp",
                  eval_grad:      "ceed_eval_grad",
                  eval_div:       "ceed_eval_div",
                  eval_curl:      "ceed_eval_curl",
                  eval_weight:    "ceed_eval_weight"}

# CeedElemTopology
shape_line     = lib.CEED_LINE
shape_triangle = lib.CEED_TRIANGLE
shape_quad     = lib.CEED_QUAD
shape_tet      = lib.CEED_TET
shape_pyramid  = lib.CEED_PYRAMID
shape_prism    = lib.CEED_PRISM
shape_hex      = lib.CEED_HEX
ceed_memtypes  = {shape_line:     "ceed_line",
                  shape_triangle: "ceed_triangle",
                  shape_quad:     "ceed_quad",
                  shape_tet:      "ceed_tet",
                  shape_pyramid:  "ceed_pyramid",
                  shape_prism:    "ceed_prism",
                  shape_hex:      "ceed_hex"}

# ------------------------------------------------------------------------------
class ceed():
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

  # Get Resource
  def getResource(self):
    """Get the full resource name for a Ceed context."""
    # libCEED call
    resource = ffi.new("char **")
    lib.CeedGetResource(self._pointer[0], resource)

    return ffi.string(resource[0]).decode("UTF-8")

  # Preferred MemType
  def getPreferredMemType(self):
    """Return Ceed preferred memory type."""
    # libCEED call
    memtype = ffi.new("CeedMemType *", ceed_mem_host)
    lib.CeedGetPreferredMemType(self._pointer[0], memtype)

    return memtype[0]

  # CeedVector
  def vector(self, size):
    """CeedVector: vector of specified length (does not allocate memory)"""
    return _Vector(self, size)

  # CeedElemRestriction

  # CeedBasis

  # CeedQFunction
  def qFunction(self, vlength, f, source):
    """CeedQFunction: independent operations at quadrature points."""
    return _QFunction(self, vlength, f, source)

  def qFunctionByName(self, name):
    """CeedQFunctionByName: independent operations at quadrature points from gallery."""
    return _QFunctionByName(self, name)

  def identityQFunction(self, size):
    """CeedIdenityQFunction: identity qfunction operation."""
    return _QFunctionIdentity(self, size)

  # CeedOperator
  def operator(self, qf, dqf, qdfT):
    """CeedOperator: composed FE-type operations on vectors."""
    return _Operator(self, qf, dqf, qdfT)

  def compositeOperator(self):
    """CompositeCeedOperator: composition of multiple CeedOperators."""
    return _CompositeOperator(self)

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedDestroy(self._pointer)

# ------------------------------------------------------------------------------
