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

  # Preferred MemType
  def getPreferredMemType(self):
    """Return Ceed preferred memory type."""
    # libCEED call
    memtype = ffi.new("CeedMemType *", lib.CEED_MEM_HOST)
    lib.CeedGetPreferredMemType(self.pointer[0], memtype)

    if (memtype[0] == lib.CEED_MEM_HOST):
      return "ceed_mem_host"
    else:
      return "ceed_mem_device"

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
