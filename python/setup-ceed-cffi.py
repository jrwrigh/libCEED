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

import os
from setuptools import setup, find_packages

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
with open(os.path.abspath("../ceed.pc.template")) as t:
  ceed_version = [line.split("Version:", 1)[1].strip() for line in t if
                    line.startswith("Version: ")]

setup(name="ceed_cffi",
      version=ceed_version[0],
      license="BSD 2",
      url="https://github.com/CEED/libCEED",
      description="libceed cffi python bindings",
      setup_requires=["cffi"],
      cffi_modules=["build-ceed-cffi.py:ffibuilder"],
      install_requires=["cffi"],
)

# ------------------------------------------------------------------------------
