# Python interface for libCEED

Preliminary development of a Python interface for libCEED. 

The file [libceed-build.py](libceed-build.py) includes the [fcci](https://cffi.readthedocs.io/en/latest/) library for the binding.

On my machine, when I run

```
python3 libceed-build.py
```

 it creates the library named `_ceed.cpython-36m-x86_64-linux-gnu.so`. I personally need to rename it to `_ceed.so` (Python doesn't like dashes in module names and I couldn't figure out how to have cffi avoid producing that suffix).

# Setup.py

Now this works and removes the middle bit

```
python setup.py build_ext
python setup.py develop

```

# Tests

The tests need to be moved up to this level to run.

# Notes

libceed.py is a shell that just pulls in the separate class files.

ceed_qfunction.py and ceed_operator.py are WIP.

# Two folder setup

I split into two packages. We build first `ceed_cffi` and then `libceed` with the commands


```
python setup.py build_ext
python setup.py develop

```
