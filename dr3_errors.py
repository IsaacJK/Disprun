#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

import os, sys

#########################################################################
# Handle Given Errors and exit program
#   Takes in a bool that tells the program to exit or not
#   Also prints out an error message if given
#########################################################################
def HandleErrors(exitBool, message):
  if exitBool == True:
    print message
    sys.exit(-1)

#########################################################################
# Iterates over data folders and checks that they all exist.
#########################################################################
def CheckFolders(parentPath, dataFolders, fitType, CheckAll=False):
  # Error handling
  errBool, missStr = False, ""
  delays = [None]

  # Check that all folders exist, otherwise just check the first
  if CheckAll == True:
    # Check that all folders exist
    missStr += "\n  Missing folders: "
    for fnum in dataFolders:
      tPath = os.path.join(parentPath, fnum)
      if not os.path.isdir(tPath):
        errBool = True
        missStr += fnum + " "

  # Check first folder for required files for peak fitting
  tPath = os.path.join(parentPath, dataFolders[0])

  # Check for SelR1rho.com file
  if not os.path.isfile(os.path.join(tPath, "SelR1rho.com")):
    missStr += "\n  Missing SelR1rho.com file"
    errBool = True

  # Check for fid.com file
  if not os.path.isfile(os.path.join(tPath, "fid.com")):
    missStr += "\n  Missing fid.com file"
    errBool = True  

  # Check for test.tab file
  if not os.path.isfile(os.path.join(tPath, "test.tab")):
    missStr += "\n  Missing test.tab file"
    errBool = True  

  # Check for vdlist
  if not os.path.isfile(os.path.join(tPath, "vdlist")):
    missStr += "\n  Missing vdlist file"
    errBool = True  
  # If vdlist exists...
  else:
    # Grab number of delays
    FILE = open(os.path.join(tPath, "vdlist"), "rU")
    delays = [x.strip() for x in FILE]
    FILE.close()

    # Check that num delays matches num ft files
    missStr += "\n  Missing ft1 files: "
    for idx, d in enumerate(delays):
      ftpath = os.path.join(tPath, "ft", "R1rho1Dpeak-%03d.ft2" % (idx+1))
      if not os.path.isfile(ftpath):
        missStr += "R1rho1Dpeak-%03d.ft2 " % (idx+1)
        errBool = True

  return errBool, missStr

