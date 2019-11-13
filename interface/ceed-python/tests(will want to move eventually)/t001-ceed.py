# @file
# Test return of CEED backend prefered memory type

import sys
import libceed

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  memtype = ceed.get_preferred_memtype()

  if (memtype == "error"):
    # LCOV_EXCL_START
    print("Error getting preferred memory type.")
    # LCOV_EXCL_STOP
