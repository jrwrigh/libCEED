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
from abc import ABC

# ------------------------------------------------------------------------------
# Ceed Enums
# ------------------------------------------------------------------------------
# CeedMemType
MEM_HOST        = lib.CEED_MEM_HOST
MEM_DEVICE      = lib.CEED_MEM_DEVICE
mem_types       = {MEM_HOST:      "ceed_mem_host",
                   MEM_DEVICE:    "ceed_mem_device"}

# CeedCopyMode
COPY_VALUES     = lib.CEED_COPY_VALUES
USE_POINTER     = lib.CEED_USE_POINTER
OWN_POINTER     = lib.CEED_OWN_POINTER
copy_modes      = {COPY_VALUES:   "ceed_copy_values",
                   USE_POINTER:   "ceed_use_pointer",
                   OWN_POINTER:   "ceed_own_pointer"}

# CeedTransposeMode
TRANSPOSE       = lib.CEED_TRANSPOSE
NOTRANSPOSE     = lib.CEED_NOTRANSPOSE
transpose_modes = {TRANSPOSE:     "ceed_transpose",
                   NOTRANSPOSE:   "ceed_notranspose"}

# CeedEvalMode
EVAL_NONE       = lib.CEED_EVAL_NONE
EVAL_INTERP     = lib.CEED_EVAL_INTERP
EVAL_GRAD       = lib.CEED_EVAL_GRAD
EVAL_DIV        = lib.CEED_EVAL_DIV
EVAL_CURL       = lib.CEED_EVAL_CURL
EVAL_WEIGHT     = lib.CEED_EVAL_WEIGHT
eval_modes      = {EVAL_NONE:     "ceed_eval_none",
                   EVAL_INTERP:   "ceed_eval_interp",
                   EVAL_GRAD:     "ceed_eval_grad",
                   EVAL_DIV:      "ceed_eval_div",
                   EVAL_CURL:     "ceed_eval_curl",
                   EVAL_WEIGHT:   "ceed_eval_weight"}

# CeedQuadMode
GAUSS           = lib.CEED_GAUSS
GAUSS_LOBATTO   = lib.CEED_GAUSS_LOBATTO
quad_modes      = {GAUSS:         "ceed_gauss",
                   GAUSS_LOBATTO: "ceed_gauss_lobatto"}

# CeedElemTopology
LINE            = lib.CEED_LINE
TRIANGLE        = lib.CEED_TRIANGLE
QUAD            = lib.CEED_QUAD
TET             = lib.CEED_TET
PYRAMID         = lib.CEED_PYRAMID
PRISM           = lib.CEED_PRISM
HEX             = lib.CEED_HEX
elem_topologies = {LINE:          "ceed_line",
                   TRIANGLE:      "ceed_triangle",
                   QUAD:          "ceed_quad",
                   TET:           "ceed_tet",
                   PYRAMID:       "ceed_pyramid",
                   PRISM:         "ceed_prism",
                   HEX:           "ceed_hex"}

# ------------------------------------------------------------------------------
# Ceed Constants
# ------------------------------------------------------------------------------

# Requests
REQUEST_IMMEDIATE = lib.CEED_REQUEST_IMMEDIATE
REQUEST_ORDERED = lib.CEED_REQUEST_ORDERED

# Object shell
class _CeedConstantObject(ABC):
  """Shell for holding constant Vector and Basis constants."""

  def __init__(self, constant):
    self._pointer = [constant]

# Vectors
VECTOR_ACTIVE = _CeedConstantObject(lib.CEED_VECTOR_ACTIVE)
VECTOR_NONE = _CeedConstantObject(lib.CEED_VECTOR_NONE)

# Basis
BASIS_COLLOCATED = _CeedConstantObject(lib.CEED_BASIS_COLLOCATED)

# ------------------------------------------------------------------------------
