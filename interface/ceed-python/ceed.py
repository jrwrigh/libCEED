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

# ------------------------------------------------------------------------------
# Ceed Enums
# ------------------------------------------------------------------------------
ceed_mem_host    = lib.CEED_MEM_HOST
ceed_mem_device  = lib.CEED_MEM_DEVICE

ceed_copy_values = lib.CEED_COPY_VALUES
ceed_use_pointer = lib.CEED_USE_POINTER
ceed_use_pointer = lib.CEED_OWN_POINTER

ceed_transpose   = lib.CEED_TRANSPOSE
ceed_notranspose = lib.CEED_NOTRANSPOSE

ceed_eval_none   = lib.CEED_EVAL_NONE
ceed_eval_interp = lib.CEED_EVAL_INTERP
ceed_eval_grad   = lib.CEED_EVAL_GRAD
ceed_eval_div    = lib.CEED_EVAL_DIV
ceed_eval_curl   = lib.CEED_EVAL_CURL
ceed_eval_weight = lib.CEED_EVAL_WEIGHT

ceed_line        = lib.CEED_LINE
ceed_triangle    = lib.CEED_TRIANGLE
ceed_quad        = lib.CEED_QUAD
ceed_tet         = lib.CEED_TET
ceed_pyramid     = lib.CEED_PYRAMID
ceed_prism       = lib.CEED_PRISM
ceed_hex         = lib.CEED_HEX

# ------------------------------------------------------------------------------
class ceed():
  """Ceed: core components."""
  # Attributes
  pointer = ffi.NULL
  vectors = []
  restrictions = []
  bases = []
  qfunctions = []
  operators = []

  # Constructor
  def __init__(self, resource = "/cpu/self"):
    # libCEED object
    self.pointer = ffi.new("Ceed *")

    # libCEED call
    resourceAscii = ffi.new("char[]", resource.encode('ascii'))
    lib.CeedInit(resourceAscii, self.pointer)

  # Get Resource
  def getResource(self):
    """Get the full resource name for a Ceed context."""
    # libCEED call
    resource = ffi.new("char **")
    lib.CeedGetResource(self.pointer[0], resource)

    return ffi.string(resource[0]).decode("UTF-8")

  # Preferred MemType
  def getPreferredMemType(self):
    """Return Ceed preferred memory type."""
    # libCEED call
    memtype = ffi.new("CeedMemType *", ceed_mem_host)
    lib.CeedGetPreferredMemType(self.pointer[0], memtype)

    return memtype[0]

  # Destructor
  def __del__(self):
    # libCEED call
    lib.CeedDestroy(self.pointer)
'''
  # CeedVector

  # CeedElemRestriction

  # CeedBasis

  # CeedQFunction
  def qFunction(self, vlength, f, source):
    # Create
    qf = QFunction(self, vlength, f, source)

    # Refrence
    self.qfunctions.append(qf)

    # Return
    return qf

  def qFunctionByName(self, name):
    # Create
    qf = QFunctionByName(self, name)

    # Refrence
    self.qfunctions.append(qf)

    # Return
    return qf

  def identityQFunction(self, size):
    # Create
    qf = QFunctionIdentity(self, size)

    # Refrence
    self.qfunctions.append(qf)

    # Return
    return qf

  # CeedOperator
  def operator(self, qf, dqf, qdfT):
    # Create
    op = Operator(self, qf, dqf, qdfT)

    # Refrence
    self.operators.append(op)

    # Return
    return op

  def compositeOperator(self):
    # Create
    op = CompositeOperator(self)

    # Refrence
    self.operators.append(op)

    # Return
    return op
'''
