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

# @file
# Test Ceed Vector functionality

import os
import libceed
import numpy as np
import check

# -------------------------------------------------------------------------------
# Utility
# -------------------------------------------------------------------------------


def check_values(ceed, x, value):
    with x.array_read() as b:
        for i in range(len(b)):
            assert b[i] == value

# -------------------------------------------------------------------------------
# Test creation, setting, reading, restoring, and destroying of a vector
# -------------------------------------------------------------------------------


def test_100(ceed_resource):
    ceed = libceed.Ceed(ceed_resource)

    n = 10
    x = ceed.Vector(n)

    a = np.arange(10, 10 + n, dtype="float64")
    x.set_array(a, cmode=libceed.USE_POINTER)

    with x.array_read() as b:
        for i in range(n):
            assert b[i] == 10 + i

# -------------------------------------------------------------------------------
# Test setValue
# -------------------------------------------------------------------------------


def test_101(ceed_resource):
    ceed = libceed.Ceed(ceed_resource)
    n = 10
    x = ceed.Vector(n)
    value = 1
    a = np.arange(10, 10 + n, dtype="float64")
    x.set_array(a, cmode=libceed.USE_POINTER)

    with x.array() as b:
        for i in range(len(b)):
            assert b[i] == 10 + i

    x.set_value(3.0)
    check_values(ceed, x, 3.0)
    del x

    x = ceed.Vector(n)
    # Set value before setting or getting the array
    x.set_value(5.0)
    check_values(ceed, x, 5.0)

# -------------------------------------------------------------------------------
# Test getArrayRead state counter
# -------------------------------------------------------------------------------


def test_102(ceed_resource):
    ceed = libceed.Ceed(ceed_resource)

    n = 10
    x = ceed.Vector(n)

    # Two read accesses should not generate an error
    a = x.get_array_read()
    b = x.get_array_read()

    x.restore_array_read()
    x.restore_array_read()

# -------------------------------------------------------------------------------
# Test setting one vector from array of another vector
# -------------------------------------------------------------------------------


def test_103(ceed_resource):
    ceed = libceed.Ceed(ceed_resource)

    n = 10

    x = ceed.Vector(n)
    y = ceed.Vector(n)

    a = np.arange(10, 10 + n, dtype="float64")
    x.set_array(a, cmode=libceed.USE_POINTER)

    with x.array() as x_array:
        y.set_array(x_array, cmode=libceed.USE_POINTER)

    with y.array_read() as y_array:
        for i in range(n):
            assert y_array[i] == 10 + i

# -------------------------------------------------------------------------------
# Test getArray to modify array
# -------------------------------------------------------------------------------


def test_104(ceed_resource):
    ceed = libceed.Ceed(ceed_resource)

    n = 10

    x = ceed.Vector(n)
    a = np.zeros(n, dtype="float64")
    x.set_array(a, cmode=libceed.USE_POINTER)

    with x.array() as b:
        b[3] = -3.14

    assert a[3] == -3.14


def test_105(ceed_resource):
    """Modification of reshaped array"""
    ceed = libceed.Ceed(ceed_resource)

    vec = ceed.Vector(12)
    with vec.array(4, 3) as x:
        x[...] = np.eye(4, 3)

    with vec.array_read(3, 4) as x:
        assert np.all(x == np.eye(4, 3).reshape(3, 4))

# -------------------------------------------------------------------------------
# Test view
# -------------------------------------------------------------------------------


def test_107(ceed_resource, capsys):
    ceed = libceed.Ceed(ceed_resource)

    n = 10
    x = ceed.Vector(n)

    a = np.arange(10, 10 + n, dtype="float64")
    x.set_array(a, cmode=libceed.USE_POINTER)

    print(x)

    stdout, stderr, ref_stdout = check.output(capsys)
    assert not stderr
    assert stdout == ref_stdout

# -------------------------------------------------------------------------------
# Test norms
# -------------------------------------------------------------------------------


def test_108(ceed_resource, capsys):
    ceed = libceed.Ceed(ceed_resource)

    n = 10
    x = ceed.Vector(n)

    a = np.arange(0, n, dtype="float64")
    for i in range(n):
        if (i % 2 == 0):
            a[i] *= -1
    x.set_array(a, cmode=libceed.USE_POINTER)

    norm = x.norm(normtype=libceed.NORM_1)

    assert abs(norm - 45.) < 1e-14

    norm = x.norm()

    assert abs(norm - np.sqrt(285.)) < 1e-14

    norm = x.norm(normtype=libceed.NORM_MAX)

    assert abs(norm - 9.) < 1e-14

# -------------------------------------------------------------------------------