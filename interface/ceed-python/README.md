## Install

To install libCEED for Python, run

```
python setup.py build
python setup.py develop --user
```

in *both* `interface/ceed-python/ceed_cffi` and `interface/ceed-python/libceed`

## Building QFunctions

To build user defined QFunctions, modify `libceed-qfunctions.c` to include
the QFunction single source file and run `python setup.qfunctions.py build`.
Renaming the output `libceed_qfunctions.*.so` may be desirable for convenience.
