import os
from cffi import FFI
ffibuilder = FFI()

Ceed = "typedef struct Ceed_private *Ceed;"

CeedInt = "typedef int32_t CeedInt;"

CeedScalar = "typedef double CeedScalar;"

CeedVector = "typedef struct CeedVector_private *CeedVector;"

CeedBasis = "typedef struct CeedBasis_private *CeedBasis;"

CeedElemRestriction = "typedef struct CeedElemRestriction_private *CeedElemRestriction;"

CeedQFunction = "typedef struct CeedQFunction_private *CeedQFunction;"

CeedOperator = "typedef struct CeedOperator_private *CeedOperator;"

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef(
  Ceed +
  """
  int CeedInit(const char *resource, Ceed *ceed);
  """
)

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("_ceed",
  """
  #include "../../include/ceed.h"   // the C header of the library
  """,
  include_dirs=[os.path.abspath('../../include')], # include path
  libraries=['ceed'],   # library name, for the linker
  library_dirs=[os.path.abspath('../../lib')], # library path, for the linker
  runtime_library_dirs=[os.path.abspath('../../lib')] # library path, at runtime
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
