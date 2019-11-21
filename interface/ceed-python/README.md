# libCEED for Python

## Install

To install libCEED for Python, run

    make

 To add instillation options, set the environment variable `PYFLAGS`, such as `PYFLAGS=--user`.

Alternatively, run

    python setup-ceed-cffi.py build_ext install
    python setup-libceed.py build install

## Building QFunctions

To build user defined QFunctions, modify `libceed-qfunctions.c` to include
the QFunction single source file and run `python setup-qfunctions.py build`.
Renaming the output `libceed_qfunctions.*.so` may be desirable for convenience.
