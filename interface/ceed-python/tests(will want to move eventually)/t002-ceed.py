# @file
# Test return of a CEED object full resource name

import sys
import libceed

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])

  resource = ceed.GetResource()

  if (resource != sys.argv[1]):
    # LCOV_EXCL_START
    print("Incorrect full resource name: " + resource + " != " + sys.argv[1])
    # LCOV_EXCL_STOP
