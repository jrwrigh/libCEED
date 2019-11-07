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

class Ceed():
  """Ceed: core components."""
  # Attributes
  self.ceed = ffi.NULL
  self.vectors = []
  self.restrictions = []
  self.bases = []
  self.qfunctions = []
  self.operators = []

  # Constructor
  def __init__(self, resource = "/cpu/self"):
    # libCEED object
    self.ceed = ffi.new("Ceed *")

    # libCEED call
    resourceAscii = resource.encode('ascii')
    libceed.CeedInit(resourceAscii, self.ceed)

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

  # Destructor
  def __def__(self):
    # libCEED call
    libceed.CeedDestroy(self.ceed)
