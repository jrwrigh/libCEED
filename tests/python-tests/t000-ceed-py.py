# @file
# Test creation and destruction of a CEED object

import sys
import libceed

if __name__ == "__main__":
  ceed = libceed.Ceed(sys.argv[1])
