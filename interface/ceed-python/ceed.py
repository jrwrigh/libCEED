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
  # Constructor
  def __init__(self, resource):
    # libCEED object
    self.ceed = ffi.new("Ceed *")

    # References to dependencies
    self.vectors = []
    self.restrictions = []
    self.bases = []
    self.qfunctions = []
    self.operators = []

    # libCEED call
    libceed.CeedInit(resource, self.ceed)

  # CeedVector

  # CeedElemRestriction

  # CeedBasis

  # CeedQFunction
  def qFunction(self, vlength, f, source):
    return QFunction(self, vlength, f, source)

  def qFunctionByName(self, name):
    return QFunctionByName(self, name)

  def identityQFunction(self, size):
    return QFunctionIdentity(self, size)

  # CeedOperator
  def operator(self, qf, dqf, qdfT):
    return Operator(self, qf, dqf, qdfT)

  def compositeOperator(self):
    return CompositeOperator(self)

  # Destructor
  def __def__(self):
    # libCEED call
    libceed.CeedDestroy(self.ceed)
