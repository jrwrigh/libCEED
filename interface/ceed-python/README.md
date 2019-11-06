# Python interface for libCEED

Preliminary development of a Python interface for libCEED. 

The file [libceed-build.py](libceed-build.py) includes the [fcci](https://cffi.readthedocs.io/en/latest/) library for the binding.

On my machine, when I run

```
python3 libceed-build.py
```

 it creates the library named `_ceed.cpython-36m-x86_64-linux-gnu.so`. I personally need to rename it to `_ceed.so` (Python doesn't like dashes in module names and I couldn't figure out how to have cffi avoid producing that suffix).
