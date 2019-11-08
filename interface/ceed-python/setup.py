from setuptools import setup, find_packages

setup(name="libceed",
      version="0.5",
      license="BSD 2",
      description="libceed python bindings",
      py_modules=["libceed"],
      setup_requires=["cffi"],
      cffi_modules=["build-ceed.py:ffibuilder"],
      install_requires=["cffi"],
)
