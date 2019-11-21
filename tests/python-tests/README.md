# libCEED Python Tests

## Building QFunctions

To build user defined QFunctions, modify `libceed-qfunctions.c` to include
the QFunction single source file and run `python setup-qfunctions.py build`.
Renaming the output `libceed_qfunctions.*.so` may be desirable for convenience.
