# @file
# Test creation and destruction of a CEED object
# \test Test creation and destruction of a CEED object

import libceed

if __name__ == "__main__":
  ceed = libceed.ceed('/cpu/self')

  memtype = ceed.getPreferredMemType()

  print(memtype)
