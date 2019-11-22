# libCEED Python Tests

## Testing

To test, first build the user QFunctions file. Then, run

  pytest test-ceed*.py --ceed /cpu/self/ref/serial

## Building QFunctions

To build user defined QFunctions, modify `libceed-qfunctions.c` to include
the QFunction single source file and run `python setup-qfunctions.py build`.
Renaming the output `libceed_qfunctions.*.so` may be desirable for convenience.
