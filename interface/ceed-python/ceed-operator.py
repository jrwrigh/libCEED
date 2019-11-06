class Operator:
    def __init__(self, ceed, qf, dqf, qdfT):
        libceed.CeedOperatorCreate(ceed, qf, dqf, dqfT, self)

    def setField(self, fieldname, restriction, lmode, basis, v):
        libceed.CeedOperatorSetField(self, fieldname, restriction,
                                     lmode, basis, v)

    def apply(self, in, out, request):
        libceed.CeedOperatorApply(self, in, out, request)

    def __del__(self):
        libceed.CeedOperatorDestroy(self)
