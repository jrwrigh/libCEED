from _ceed import ffi, lib

if __name__ == "__main__":
  ceed = ffi.new("Ceed *")
  resource = ffi.new("char[]", "/cpu/self".encode('ascii'))
  lib.CeedInit(resource, ceed)
