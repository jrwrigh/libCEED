from _ceed import ffi, lib

if __name__ == "__main__":
  ceed = ffi.new("Ceed")
  lib.CeedInit("cpu/self/ref/serial", ceed)
  print(ceed)
  print(ceed[0])
