# libCEED for Python

## Install

To install libCEED for Python, run

    make

To add installation options, set the environment variable `PYFLAGS`, such as `PYFLAGS=--user`.

Alternatively, run

    python setup-ceed-cffi.py build_ext install
    python setup-libceed.py build install
