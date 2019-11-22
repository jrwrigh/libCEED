# @file
# Test grad with a 2D Simplex non-tensor H1 basis

import sys, math
import libceed
import numpy as np
import buildmats as bm

def feval(x1, x2):
  return x1*x1 + x2*x2 + x1*x2 + 1

def dfeval(x1, x2):
  return 2*x1 + x2

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  P, Q, dim = 6, 4, 2

  xq = np.array([0.2, 0.6, 1./3., 0.2, 0.2, 0.2, 1./3., 0.6], dtype="float64")
  xr = np.array([0., 0.5, 1., 0., 0.5, 0., 0., 0., 0., 0.5, 0.5, 1.], dtype="float64")
  in_array = np.empty(P, dtype="float64")
  qref = np.empty(dim*Q, dtype="float64")
  qweight = np.empty(Q, dtype="float64")

  interp, grad = bm.buildmats(qref, qweight)

  b = ceed.BasisH1(libceed.TRIANGLE, 1, P, Q, interp, grad, qref, qweight)

  # Interpolate function to quadrature points
  for i in range(P):
    in_array[i] = feval(xr[0*P+i], xr[1*P+i])

  in_vec = ceed.Vector(P)
  in_vec.set_array(in_array, cmode=libceed.USE_POINTER)
  out_vec = ceed.Vector(Q*dim)
  out_vec.set_value(0)

  b.apply(1, libceed.EVAL_GRAD, in_vec, out_vec)

  # Check values at quadrature points
  out_array = out_vec.get_array_read()
  for i in range(Q):
    value = dfeval(xq[0*Q+i], xq[1*Q+i])
    if math.fabs(out_array[0*Q+i] - value) > 1E-10:
      # LCOV_EXCL_START
      print("[%d] %f != %f"%(i, out[0*Q+i], value))
    # LCOV_EXCL_STOP
    value = dfeval(xq[1*Q+i], xq[0*Q+i])
    if math.fabs(out_array[1*Q+i] - value) > 1E-10:
      # LCOV_EXCL_START
      print("[%d] %f != %f"%(i, out[1*Q+i], value))
    # LCOV_EXCL_STOP

  out_vec.restore_array_read()
