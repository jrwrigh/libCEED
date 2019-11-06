class QFunction:
    def __init__(self, ceed, ):

    def addInput(self, fieldname, size, emode):
        libceed.CeedQFunctionAddInput(self, fieldname, size, emode)

    def addOutput(self, fieldname, size, emode):
        libceed.CeedQFunctionAddOutput(self, fieldname, size, emode)

    def apply(self, q, u, v):
        libceed.CeedQFunctionApply(self, q, u, v)

    def __del__(self):
        libceed.CeedQFunctionDestroy(self)
