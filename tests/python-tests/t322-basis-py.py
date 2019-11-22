# @file
# Test integration with a 2D Simplex non-tensor H1 basis

import sys, math
import libceed
import numpy as np
import buildmats as bm

def feval(x1, x2):
  return x1*x1 + x2*x2 + x1*x2 + 1

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  P, Q, dim = 6, 4, 2

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
  out_vec = ceed.Vector(Q)
  out_vec.set_value(0)
  weights_vec = ceed.Vector(Q)
  weights_vec.set_value(0)

  b.apply(1, libceed.EVAL_INTERP, in_vec, out_vec)
  b.apply(1, libceed.EVAL_WEIGHT, libceed.VECTOR_NONE, weights_vec)

  # Check values at quadrature points
  out_array = out_vec.get_array_read()
  weights_array = weights_vec.get_array_read()
  sum = 0
  for i in range(Q):
    sum += out_array[i]*weights_array[i]
  if math.fabs(sum - 17./24.) > 1E-10:
    print("%f != %f"%(sum, 17./24.))

  out_vec.restore_array_read()
  weights_vec.restore_array_read()
