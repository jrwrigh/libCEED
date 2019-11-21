# libCEED for Python

## Install

To install libCEED for Python, run

    make python

in the `libCEED` directory, or run

```
python setup.py build
python setup.py install
```

in *both* `interface/ceed-python/ceed-cffi` and `interface/ceed-python/libceed`

## Building QFunctions

To build user defined QFunctions, modify `libceed-qfunctions.c` to include
the QFunction single source file and run `python setup-qfunctions.py build`.
Renaming the output `libceed_qfunctions.*.so` may be desirable for convenience.
