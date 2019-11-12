# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

import os
from cffi import FFI
ffibuilder = FFI()

# ------------------------------------------------------------------------------
# Strings for cdef
# ------------------------------------------------------------------------------
CeedInt = "typedef int32_t CeedInt;"

CeedScalar = "typedef double CeedScalar;"

Ceed = "typedef struct Ceed_private *Ceed;"

CeedRequest = "typedef struct CeedRequest_private *CeedRequest;"

CeedVector = "typedef struct CeedVector_private *CeedVector;"

CeedBasis = "typedef struct CeedBasis_private *CeedBasis;"

CeedElemRestriction = "typedef struct CeedElemRestriction_private *CeedElemRestriction;"

CeedQFunction = "typedef struct CeedQFunction_private *CeedQFunction;"

CeedOperator = "typedef struct CeedOperator_private *CeedOperator;"

CeedMemType = "typedef enum { CEED_MEM_HOST, CEED_MEM_DEVICE, } CeedMemType;"

CeedCopyMode = "typedef enum { CEED_COPY_VALUES, CEED_USE_POINTER, CEED_OWN_POINTER, } CeedCopyMode;"

CeedTransposeMode = "typedef enum { CEED_NOTRANSPOSE, CEED_TRANSPOSE } CeedTransposeMode;"

CeedEvalMode = "typedef enum { CEED_EVAL_NONE   = 0, CEED_EVAL_INTERP = 1, CEED_EVAL_GRAD   = 2, CEED_EVAL_DIV    = 4, CEED_EVAL_CURL   = 8, CEED_EVAL_WEIGHT = 16} CeedEvalMode;"

CeedQuadMode = "typedef enum { CEED_GAUSS = 0, CEED_GAUSS_LOBATTO = 1, } CeedQuadMode;"

CeedElemTopology = "typedef enum { CEED_LINE = 1 << 16 | 0, CEED_TRIANGLE = 2 << 16 | 1, CEED_QUAD = 2 << 16 | 2, CEED_TET = 3 << 16 | 3, CEED_PYRAMID = 3 << 16 | 4, CEED_PRISM = 3 << 16 | 5, CEED_HEX = 3 << 16 | 6, } CeedElemTopology;"

# ------------------------------------------------------------------------------
# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
# ------------------------------------------------------------------------------
ffibuilder.cdef(
  CeedInt +
  CeedScalar +
  CeedMemType +
  CeedCopyMode +
  CeedTransposeMode +
  CeedEvalMode +
  CeedQuadMode +
  CeedElemTopology +
  Ceed +
  """
  int CeedInit(const char *resource, Ceed *ceed);
  int CeedGetResource(Ceed ceed, const char **resource);
  int CeedDestroy(Ceed *ceed);
  int CeedGetPreferredMemType(Ceed ceed, CeedMemType *type);
  """ +
  CeedRequest +
  """
  CeedRequest *const CEED_REQUEST_IMMEDIATE;
  CeedRequest *const CEED_REQUEST_ORDERED;
  int CeedRequestWait(CeedRequest *req);
  """ +
  CeedVector +
  """
  extern CeedVector CEED_VECTOR_ACTIVE;
  extern CeedVector CEED_VECTOR_NONE;
  int CeedVectorCreate(Ceed ceed, CeedInt len, CeedVector *vec);
  int CeedVectorSetArray(CeedVector vec, CeedMemType mtype,
                                   CeedCopyMode cmode, CeedScalar *array);
  int CeedVectorSetValue(CeedVector vec, CeedScalar value);
  int CeedVectorSyncArray(CeedVector vec, CeedMemType mtype);
  int CeedVectorGetArray(CeedVector vec, CeedMemType mtype,
                                   CeedScalar **array);
  int CeedVectorGetArrayRead(CeedVector vec, CeedMemType mtype,
                                       const CeedScalar **array);
  int CeedVectorRestoreArray(CeedVector vec, CeedScalar **array);
  int CeedVectorRestoreArrayRead(CeedVector vec,
    const CeedScalar **array);
  int CeedVectorView(CeedVector vec, const char *fpfmt, FILE *stream);
  int CeedVectorGetLength(CeedVector vec, CeedInt *length);
  int CeedVectorDestroy(CeedVector *vec);
  """ +
  CeedElemRestriction +
  """
  int CeedElemRestrictionCreate(Ceed ceed, CeedInt nelem,
    CeedInt elemsize, CeedInt nnodes, CeedInt ncomp, CeedMemType mtype,
    CeedCopyMode cmode,
    const CeedInt *indices, CeedElemRestriction *rstr);
  int CeedElemRestrictionCreateIdentity(Ceed ceed, CeedInt nelem,
    CeedInt elemsize, CeedInt nnodes, CeedInt ncomp, CeedElemRestriction *rstr);
  int CeedElemRestrictionCreateBlocked(Ceed ceed, CeedInt nelem,
    CeedInt elemsize, CeedInt blksize, CeedInt nnodes, CeedInt ncomp,
    CeedMemType mtype,
    CeedCopyMode cmode, const CeedInt *indices, CeedElemRestriction *rstr);
  int CeedElemRestrictionCreateVector(CeedElemRestriction rstr,
    CeedVector *lvec, CeedVector *evec);
  int CeedElemRestrictionApply(CeedElemRestriction rstr,
    CeedTransposeMode tmode, CeedTransposeMode lmode, CeedVector u,
    CeedVector ru, CeedRequest *request);
  int CeedElemRestrictionApplyBlock(CeedElemRestriction rstr,
    CeedInt block, CeedTransposeMode tmode, CeedTransposeMode lmode,
    CeedVector u, CeedVector ru, CeedRequest *request);
  int CeedElemRestrictionGetMultiplicity(CeedElemRestriction rstr,
    CeedVector mult);
  int CeedElemRestrictionView(CeedElemRestriction rstr, FILE *stream);
  int CeedElemRestrictionDestroy(CeedElemRestriction *rstr);
  """ +
  CeedBasis +
  """
  extern CeedBasis CEED_BASIS_COLLOCATED;
  int CeedBasisCreateTensorH1Lagrange(Ceed ceed, CeedInt dim,
    CeedInt ncomp, CeedInt P, CeedInt Q, CeedQuadMode qmode, CeedBasis *basis);
  int CeedBasisCreateTensorH1(Ceed ceed, CeedInt dim, CeedInt ncomp,
                                        CeedInt P1d, CeedInt Q1d,
                                        const CeedScalar *interp1d,
                                        const CeedScalar *grad1d,
                                        const CeedScalar *qref1d,
                                        const CeedScalar *qweight1d,
                                        CeedBasis *basis);
  int CeedBasisCreateH1(Ceed ceed, CeedElemTopology topo,
                                  CeedInt ncomp,
                                  CeedInt nnodes, CeedInt nqpts,
                                  const CeedScalar *interp,
                                  const CeedScalar *grad,
                                  const CeedScalar *qref,
                                  const CeedScalar *qweight, CeedBasis *basis);
  int CeedBasisView(CeedBasis basis, FILE *stream);
  int CeedBasisGetNumNodes(CeedBasis basis, CeedInt *P);
  int CeedBasisGetNumQuadraturePoints(CeedBasis basis, CeedInt *Q);
  int CeedBasisApply(CeedBasis basis, CeedInt nelem,
                               CeedTransposeMode tmode,
                               CeedEvalMode emode, CeedVector u, CeedVector v);
  int CeedBasisDestroy(CeedBasis *basis);

  int CeedGaussQuadrature(CeedInt Q, CeedScalar *qref1d,
                                    CeedScalar *qweight1d);
  int CeedLobattoQuadrature(CeedInt Q, CeedScalar *qref1d,
                                      CeedScalar *qweight1d);
  int CeedQRFactorization(Ceed ceed, CeedScalar *mat, CeedScalar *tau,
                                    CeedInt m, CeedInt n);
  int CeedSymmetricSchurDecomposition(Ceed ceed, CeedScalar *mat,
    CeedScalar *lambda, CeedInt n);
  int CeedSimultaneousDiagonalization(Ceed ceed, CeedScalar *matA,
    CeedScalar *matB, CeedScalar *x, CeedScalar *lambda, CeedInt n);
  """ +
  CeedQFunction +
  """

typedef int (*CeedQFunctionUser)(void *ctx, const CeedInt Q,
                                 const CeedScalar *const *in,
                                 CeedScalar *const *out);

  int CeedQFunctionCreateInterior(Ceed ceed, CeedInt vlength,
    CeedQFunctionUser f, const char *source, CeedQFunction *qf);
  int CeedQFunctionCreateInteriorByName(Ceed ceed, const char *name,
    CeedQFunction *qf);
  int CeedQFunctionCreateIdentity(Ceed ceed, CeedInt size,
    CeedQFunction *qf);
  int CeedQFunctionAddInput(CeedQFunction qf, const char *fieldname,
                                      CeedInt size, CeedEvalMode emode);
  int CeedQFunctionAddOutput(CeedQFunction qf, const char *fieldname,
                                       CeedInt size, CeedEvalMode emode);
  int CeedQFunctionSetContext(CeedQFunction qf, void *ctx,
                                        size_t ctxsize);
  int CeedQFunctionApply(CeedQFunction qf, CeedInt Q,
                                   CeedVector *u, CeedVector *v);
  int CeedQFunctionDestroy(CeedQFunction *qf);
  """ +
  CeedOperator +
  """
  int CeedOperatorCreate(Ceed ceed, CeedQFunction qf,
                                   CeedQFunction dqf, CeedQFunction dqfT,
                                   CeedOperator *op);
  int CeedCompositeOperatorCreate(Ceed ceed, CeedOperator *op);
  int CeedOperatorSetField(CeedOperator op, const char *fieldname,
                                     CeedElemRestriction r,
                                     CeedTransposeMode lmode, CeedBasis b,
                                     CeedVector v);
  int CeedCompositeOperatorAddSub(CeedOperator compositeop,
    CeedOperator subop);
  int CeedOperatorAssembleLinearQFunction(CeedOperator op,
    CeedVector *assembled, CeedElemRestriction *rstr, CeedRequest *request);
  int CeedOperatorAssembleLinearDiagonal(CeedOperator op,
    CeedVector *assembled, CeedRequest *request);
  int CeedOperatorApply(CeedOperator op, CeedVector in,
                                  CeedVector out, CeedRequest *request);
  int CeedOperatorDestroy(CeedOperator *op);
  """
)

# ------------------------------------------------------------------------------
# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
# ------------------------------------------------------------------------------
ffibuilder.set_source("_ceed",
  """
  #include <ceed.h>   // the C header of the library
  """,
  include_dirs=[os.path.abspath("../../include")], # include path
  libraries=["ceed"],   # library name, for the linker
  library_dirs=[os.path.abspath("../../lib")], # library path, for the linker
  runtime_library_dirs=[os.path.abspath("../../lib")] # library path, at runtime
)

# ------------------------------------------------------------------------------
# Builder
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

# ------------------------------------------------------------------------------
