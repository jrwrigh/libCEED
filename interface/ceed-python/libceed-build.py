import os
from cffi import FFI
ffibuilder = FFI()

objdelegate = """typedef struct {
  char *objname;
  Ceed delegate;
} objdelegate;"""

CeedMemType = """typedef enum {
  /// Memory resides on the host
  CEED_MEM_HOST,
  /// Memory resides on a device (corresponding to \ref Ceed resource)
  CEED_MEM_DEVICE,
} CeedMemType;"""

CeedInt = "typedef int32_t CeedInt;"

Ceed_private = """struct Ceed_private {
  const char *resource;
  Ceed delegate;
  Ceed parent;
  objdelegate *objdelegates;
  int objdelegatecount;
  int (*Error)(Ceed, const char *, int, const char *, int, const char *,
               va_list);
  int (*GetPreferredMemType)(CeedMemType *);
  int (*Destroy)(Ceed);
  int (*VectorCreate)(CeedInt, CeedVector);
  int (*ElemRestrictionCreate)(CeedMemType, CeedCopyMode,
                               const CeedInt *, CeedElemRestriction);
  int (*ElemRestrictionCreateBlocked)(CeedMemType, CeedCopyMode,
                                      const CeedInt *, CeedElemRestriction);
  int (*BasisCreateTensorH1)(CeedInt, CeedInt, CeedInt, const CeedScalar *,
                             const CeedScalar *, const CeedScalar *,
                             const CeedScalar *, CeedBasis);
  int (*BasisCreateH1)(CeedElemTopology, CeedInt, CeedInt, CeedInt,
                       const CeedScalar *,
                       const CeedScalar *, const CeedScalar *,
                       const CeedScalar *, CeedBasis);
  int (*TensorContractCreate)(CeedBasis, CeedTensorContract);
  int (*QFunctionCreate)(CeedQFunction);
  int (*OperatorCreate)(CeedOperator);
  int (*CompositeOperatorCreate)(CeedOperator);
  int refcount;
  void *data;
  foffset *foffsets;
};"""

Ceed = "typedef struct Ceed_private *Ceed;"

CeedVector = "typedef struct CeedVector_private *CeedVector;"

CeedCopyMode = """typedef enum {
  CEED_COPY_VALUES,
  CEED_USE_POINTER,
  CEED_OWN_POINTER,
} CeedCopyMode;"""

CeedElemRestriction = "typedef struct CeedElemRestriction_private *CeedElemRestriction;"

CeedScalar = "typedef double CeedScalar;"

CeedBasis = "typedef struct CeedBasis_private *CeedBasis;"

CeedElemTopology = """typedef enum {
  CEED_LINE = 1 << 16 | 0,
  CEED_TRIANGLE = 2 << 16 | 1,
  CEED_QUAD = 2 << 16 | 2,
  CEED_TET = 3 << 16 | 3,
  CEED_PYRAMID = 3 << 16 | 4,
  CEED_PRISM = 3 << 16 | 5,
  CEED_HEX = 3 << 16 | 6,
} CeedElemTopology;"""

CeedTensorContract = "typedef struct CeedTensorContract_private *CeedTensorContract;"

CeedQFunction = "typedef struct CeedQFunction_private *CeedQFunction;"

CeedOperator = "typedef struct CeedOperator_private *CeedOperator;"

foffset = """typedef struct {
  const char *fname;
  size_t offset;
} foffset;"""

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef(
  "typedef ... va_list;" +
  Ceed +
  objdelegate +
  CeedMemType +
  CeedInt +
  CeedVector +
  CeedCopyMode +
  CeedElemRestriction +
  CeedScalar +
  CeedBasis +
  CeedElemTopology +
  CeedTensorContract +
  CeedQFunction +
  CeedOperator +
  foffset +
  Ceed_private +
  """
  int CeedInit(const char *resource, Ceed *ceed);
""")

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("_ceed",
"""
  #include "../../include/ceed.h"   // the C header of the library
  #include "../../include/ceed-impl.h"   // the C header of the library
""",
  include_dirs=[os.path.abspath('../../include')], # include path
  libraries=['ceed'],   # library name, for the linker
  library_dirs=[os.path.abspath('../../lib')], # library path, for the linker
  runtime_library_dirs=[os.path.abspath('../../lib')] # library path, at runtime
  )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
