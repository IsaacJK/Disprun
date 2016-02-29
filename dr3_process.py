#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

import os, sys
import nmrglue as ng
from numpy import log10, array, insert, sort

def makeFolder(pathToFolder):
  if not os.path.exists(pathToFolder): 
    os.makedirs(pathToFolder) 
#---------------------------#---------------------------#
# Goes in to noisePath and reads out noise
#  corresponding to delays in delay list
# Returns array with [[delay, noise]...]
def GetNoise(noisePath, delays, fitType):
  # Noise list
  noiseList = []
  
  if fitType == "1d":
    # Loop over delays to read out noise
    for idx, i in enumerate(delays):
      noiseFile = "R1rho1Dnoise-%03d.txt" % (idx+1)
      tNpath = os.path.join(noisePath, noiseFile)
      FILE = open(tNpath, "rU")
      noise = [i] + [float(x.strip()) for x in FILE]
      FILE.close()
      noiseList.append(noise)

  elif fitType == "2d":
    # Loop over delays to read out noise
    for idx, i in enumerate(delays):
      noiseFile = "R1rho1Dnoise-001.txt"
      tNpath = os.path.join(noisePath, noiseFile)
      FILE = open(tNpath, "rU")
      noise = [i] + [float(x.strip()) for x in FILE]
      FILE.close()
      noiseList.append(noise) 

  return array(noiseList)

#---------------------------#---------------------------#
# Reads in a list of folders, grabs R1rho relevant
#  parameters from acqus file and returns a dict of lists
#  containing these parameters
# f.num : 0 cnst28, 1 inv.cnst28, 2 cnst12, 3 pl23,
#   4 te, 5 o1p, 6 o2p, 7 o3p, 8 p1, 9 nbl, 10 d30,
#   11 d31, 12 vdlist, 13 pulprog
#---------------------------#---------------------------#
def GrabSLPOffs(acqusPath, folderNum):
  # Data to return
  slpoffs = {}

  # Read acqus files
  tdict, tdata = ng.fileio.bruker.read(dir=acqusPath, read_pulseprogram=False)
  # Base freqs
  bf1 = tdict['acqus']['BF1']
  bf2 = tdict['acqus']['BF2']
  bf3 = tdict['acqus']['BF3']
  # w1h
  cnst11 = tdict['acqus']['CNST'][11]
  # SLP
  cnst12 = tdict['acqus']['CNST'][12]
  cnst28 = tdict['acqus']['CNST'][28]
  cnst28i = -1*cnst28
  cnst30 = tdict['acqus']['CNST'][30]
  o1p = tdict['acqus']['O1'] / tdict['acqus']['BF1']
  o2p = tdict['acqus']['O2'] / tdict['acqus']['BF2']
  o3p = tdict['acqus']['O3'] / tdict['acqus']['BF3']
  ns = tdict['acqus']['NS']
  ds = tdict['acqus']['DS']
  rg = tdict['acqus']['RG']
  ds = tdict['acqus']['DS']
  nbl = tdict['acqus']['NBL']
  # 1H hard 90
  p1 = tdict['acqus']['P'][1]
  # 13C hard 90
  p3 = tdict['acqus']['P'][3]
  # 15N hard 90
  p21 = tdict['acqus']['P'][21]
  # 15N sinc pulse on water
  p25 = tdict['acqus']['P'][25]
  # Decoupling pulses
  # for carbon cpd dec during Trelax
  pcpd1 = tdict['acqus']['PCPD'][1]
  # 13C dec
  pcpd2 = tdict['acqus']['PCPD'][2]
  # 15N dec
  pcpd3 = tdict['acqus']['PCPD'][3]
  d30 = tdict['acqus']['D'][30]
  d31 = tdict['acqus']['D'][31]
  te = tdict['acqus']['TE']
  vdlist = tdict['acqus']['VDLIST']
  pp = tdict['acqus']['PULPROG']
  sp5 = tdict['acqus']['SP'][5]
  try:
    pl12 = tdict['acqus']['PLW'][12]
    pl16 = tdict['acqus']['PLW'][16]
    pl23 = tdict['acqus']['PLW'][23]
    pl24 = tdict['acqus']['PLW'][24]
    # Convert to att. dB
    pl12 = round(-10.*log10(pl12),2)
    pl16 = round(-10.*log10(pl16),2)
    pl23 = round(-10.*log10(pl23),2)
    pl24 = round(-10.*log10(pl24),2)
  except KeyError:
    pl12 = tdict['acqus']['PL'][12]
    pl16 = tdict['acqus']['PL'][16]
    pl23 = tdict['acqus']['PL'][23]
    pl24 = tdict['acqus']['PL'][24]    
  slpoffs = {"FN": folderNum, "cnst28": cnst28, 
             "cnst28i": cnst28i, "cnst12": cnst12,
             "pl23": pl23, "cnst11": cnst11, "pl24": pl24, 
             "te": te, "o1p": o1p, "o2p": o2p, "o3p": o3p, "cnst30": cnst30,
             "p1": p1, "p3": p3, "p21": p21,
             "p25": p25, "sp5": sp5, "ns": ns, "ds": ds, "nbl": nbl, "d30": d30, "d31": d31,
             "vdlist": vdlist, "rg": rg, "pcpd1": pcpd1, "pcpd2": pcpd2, "pldb12": pl12,
             "pcpd3": pcpd3, "pldb16": pl16, "pulprog": pp, "bf1":bf1, "bf2": bf2, "bf3": bf3}
  # slpoffs = {"Folder Number": folderNum, "Offset (Hz cnst28)": cnst28, 
  #            "Corr. Offset (Hz cnst28 inv)": cnst28i, "SLP (Hz cnst12)": cnst12,
  #            "SLP (dB pl23)": pl23, "w1H (Hz cnst11)": cnst11, "w1H (dB pl24)": pl24, 
  #            "Temp (K)": te, "O1P (ppm)": o1p, "O2P (ppm)": o2p, "O3P (ppm)": o3p,
  #            "P1 1H 90o (us)": p1, "P3 13C 90o (us)": p3, "P21 15N 90o (us)": p21,
  #            "P25 (us)": p25, "SP5": sp5, "NS": ns, "DS": ds, "NBL": nbl, "D30": d30, "D31": d31,
  #            "VDlist": vdlist, "RG": rg, "PCPD1": pcpd1, "PCPD2": pcpd2, "PL12 (dB)": pl12,
  #            "PCPD3": pcpd3, "PL16 (dB)": pl16, "pulprog": pp}

  return slpoffs

#---------------------------#---------------------------#
# Goes in to specified path and removes all 'ms' or
#  spaces in vdlist
#---------------------------#---------------------------#
def CleanDly(ppath):
  vdlistpath = os.path.join(ppath, "vdlist")
  if os.path.isfile(vdlistpath):
    # Strip spaces from vdlist
    FILE = open(vdlistpath, "rU")
    delays = [x.strip() for x in FILE if len(x.strip()) != 0]
    FILE.close()
    # check to see if ms in vdlist
    if "ms" in delays[0]:
      delays = [str(float(x.replace("ms", ""))/1e3) for x in delays if len(x.replace("ms", "")) > 0]
    FILE = open(vdlistpath, "wb")
    # Write out cleaned vdlist
    for line in delays:
      FILE.write(line + "\n")
    FILE.close()

#---------------------------#---------------------------#
# Reads updated delay list (tau.list), noise,
#  and intensities from outT.tab
# Returns 2 numpy arrays of delays, intensities and noise
#  1. raw intensities and error
#  2. normalized intensities and error
#---------------------------#---------------------------#
def GrabDlyIntsNoise(ParentPath, peakNum, fitType):
  # List for noise
  noiseList = []
  # List for intensities
  intsList = []

  # Set delay path
  # if fitType == "1d":
    # Get delays from vdlist for 1d
  tauPath = os.path.join(ParentPath, "vdlist")
  # elif fitType == "2d":
  #   # Get delays from vdlist for 1d
  #   tauPath = os.path.join(ParentPath, "tau.list")
  FILE = open(tauPath, "rU")
  newDelays = array([float(x.strip()) for x in FILE if len(x.strip()) > 0])
  FILE.close()

  
  if fitType == "2d":
    newDelays = array(sorted(list(newDelays)))
  # Legnth of delay list
  numdly = len(newDelays)

  # Get the noise value from the R1rho1Dnoise*.txt files
  #  Noise for 1D comes from RMS for each 1D slice
  #  Noise for 2D is just a list of the same noise
  #   from the first slice
  noiseList = GetNoise(ParentPath, newDelays, fitType)

  # Grab intensities for peaks
  if fitType == "1d":
    for idx in range(len(newDelays)):
      # Grab intensities correlating to these delays
      tabpath = os.path.join(ParentPath, "R1rho1DFit-%03d.tab" % (idx+1)) 
      pc,pf,rec = ng.fileio.pipe.read_table(tabpath)
      for line in rec:
        if str(line[0]) == peakNum:
          # Cast line as list
          line = list(line)

          # Read last N columns of list, these are intensities
          ints = array(line[-numdly:])

          # Grab noise, and normalize it to raw intensity
          rawInt = float(line[9])
          intsList.append(rawInt) 

  elif fitType == "2d":
    # Grab intensities correlating to these delays
    tabpath = os.path.join(ParentPath, "R1rho1DFit-001.tab")
    pc,pf,rec = ng.fileio.pipe.read_table(tabpath)
    for line in rec:
      if str(line[0]) == peakNum:
        # Cast line as list
        line = list(line)

        # Read last N columns of list, these are intensities
        intsList = array(line[-numdly:])

        # Grab noise, and normalize it to raw intensity
        rawInt = float(line[9])

        # Calculate relative noise
        noiseList = array([array([x,ye/rawInt])
                           for x,ye in
                           zip(noiseList[:,0], noiseList[:,1])])
  
  ### Arrange [delay, intensity, error] numpy arrays
  # rDIN = [delay, raw intensity, raw error]
  # nDIN = [delay, normalized intensity, normalized error]
  # For pseudo-2D version, rDIN == nDIN
  # Insert intensities into delay numpy array
  # [dly, intensity, noise]
  intsList = array(intsList)
  rDIN = insert(noiseList,1,intsList,axis=1)
  # Normalize intensities to 1, then normalize noise to intensity
  nDIN = array([array([x, y/rDIN[:,1].max(), ye/y])
                    for x,y,ye in zip(rDIN[:,0],rDIN[:,1],rDIN[:,2])])
  # Sort the lists
  rDIN = rDIN[rDIN[:,0].argsort()]
  nDIN = nDIN[nDIN[:,0].argsort()]

  return rDIN, nDIN

#---------------------------#---------------------------#
# Generates a SelR1rho.com script in path defined
#  Returns the path to that script
# Adapted from Yi's 'get1d.com' script
#---------------------------#---------------------------#
def WriteSelCom(WritePath, delays, fitType, p0=0.0):
  # Output path and name for com script
  outSelR1rho = os.path.join(WritePath, "SelR1rho.com")
  outFit = os.path.join(WritePath, "fit.com")
  dlyStr = ""
  for d in delays:
    dlyStr += str(d) + " "

  if fitType == "2d":
    # Construct Comscript
    SelCom = '''
#!/bin/csh

set tauList = (%s)

nmrPipe -in test.fid            \\
| nmrPipe  -fn SOL            \\
| nmrPipe  -fn POLY -time           \\
| nmrPipe  -fn SP -off 0.5 -end 0.95 -pow 4 -c 0.5    \\
| nmrPipe  -fn ZF -auto -zf 4         \\
| nmrPipe  -fn FT -auto           \\
| nmrPipe  -fn PS -p0 %s -p1 0.00  -verb     \\
| nmrPipe  -fn POLY -auto -ord 4        \\
   -ov -out out.ft2
''' % (dlyStr, p0)
    SelCom += "split2D.com -in out.ft2 -outDir ft -outName R1rho1Dpeak-%03d.ft2 -tau $tauList\n"

    SelCom += '''
if (-e tau.list) then
   /bin/rm tau.list
endif

foreach f (ft/R1rho1Dpeak-*.ft2)
   set tau = (`getParm -in $f -parm FDTAU`)
   echo $tau >> tau.list
end
'''

    FILE = open(outSelR1rho, "wb")
    FILE.writelines(SelCom)
    FILE.close()

    FitCom = '''
#!/bin/csh

autoFit.tcl -specName ft/R1rho1Dpeak-%03d.ft2 -inTab test.tab -series \\
            -modX LORENTZ1D -outTab R1rho1DFit-001.tab
awk '/^set noiseRMS/ {print $4}' autoFit.com > R1rho1Dnoise-001.txt
rm autoFit.com
rm axt.tab
'''
  
    FILE = open(outFit, "wb")
    FILE.writelines(FitCom)
    FILE.close()

  else:
    # Make ft folder
    makeFolder(os.path.join(WritePath, "ft"))
    SelCom = '''
#!/bin/csh

set tauList = (%s)

nmrPipe -in test.fid            \\
| nmrPipe  -fn SOL            \\
| nmrPipe  -fn POLY -time           \\
| nmrPipe  -fn SP -off 0.5 -end 0.95 -pow 3 -c 0.5    \\
| nmrPipe  -fn ZF -auto -zf 4         \\
| nmrPipe  -fn FT -auto           \\
| nmrPipe  -fn PS -p0 %s -p1 0.00  -verb     \\
| nmrPipe  -fn POLY -auto -ord 4        \\
   -ov -out out.ft2
''' % (dlyStr, p0)
  
    for idx in range(len(delays)):
      SelCom += "readROI -in out.ft2 -ndim 1 -x 1H -dy F1 %03d -out ft/R1rho1Dpeak-%03d.ft2\n" % (idx+1, idx+1)

    FILE = open(outSelR1rho, "wb")
    FILE.writelines(SelCom)
    FILE.close()

    FitCom = '''
#!/bin/csh

rm -f aux.tab
'''
    for idx in range(len(delays)):
      FitCom += '''
autoFit.tcl -specName ft/R1rho1Dpeak-%03d.ft2 -inTab test.tab -outTab R1rho1DFit-%03d.tab
awk '/^set noiseRMS/ {print $4}' autoFit.com > R1rho1Dnoise-%03d.txt
rm autoFit.com
rm axt.tab
''' % (idx+1, idx+1, idx+1)
  
    FILE = open(outFit, "wb")
    FILE.writelines(FitCom)
    FILE.close()

  return outSelR1rho


#---------------------------#---------------------------#
# Takes a path to an input text file containing:
#   PeakNum [Fit peak number, e.g. 1]
#   DataFolders [Folder range, e.g. 1400-1539]
#   OutputFolder [Folder to dump out fits, e.g. "R1rho"]
# Return error flags, and
#  peakNum : int to mark which peak to choose for fitting
#  dataFolders : numpy array of N data folders (numerical)
#  outPath : Name of directory to place fits
#  fitType : '1d' or '2d', for 1d or pseudo-2d
#  errType : 'std', 'mc', or 'bs' for std err, mc-err, bootstrap-err
#---------------------------#---------------------------#
def ParseInp(inpPath):
  # Error flags
  errBool, missStr = False, ""
  # To store retrieved values
  peakNum, dataFolders, outPath = None, None, None
  fitType, readFlag = "1d", "no"
  # Variables in input file
  variables = set(["peaknum", "datafolders", "outputfolder", "fittype", "read"])

  FILE = open(inpPath, "rU")
  # rawdata = [x.strip().split() for x in FILE]
  rawdata = [x.strip().split() for x in FILE 
             if len(x.strip().split()) > 0 and
             "#" not in x.strip().split() and
             x.strip().split()[0].lower() in variables]
  FILE.close()

  # Check that len is == 3, else, flag error
  if len(rawdata) < 4:
    missStr += "\n  Too few parameters in input text file."
    errBool = True
  elif len(rawdata) > 4:
    missStr += "\n  Too many parameters in input text file."
    errBool = True   
  
  for val in rawdata:
    # Check peak number, assign if correct
    if val[0].lower() == "peaknum":
      try:
        # Needs to be numeric
        int(val[1])
        peakNum = val[1]
      except ValueError:
        missStr += "\n  Peak number is non-numeric (%s)" % val[1]
        errBool = True
    # Convert range of numbers to a list of numbers,
    #   i.e. 1-3 to ['1', '2', '3']
    elif val[0].lower() == "datafolders":
      # Return sorted list of numbers, and error flags
      errBool, tstr, dataFolders = ParseIndexString(val[1])
      # Update miss-str
      missStr += tstr
    # No checks here
    elif val[0].lower() == "outputfolder":
      outPath = val[1]
    # Check for peak fit type flag
    elif val[0].lower() == "fittype":
      if "1d" in val[1].lower():
        fitType = "1d"
      elif "2d" in val[1].lower():
        fitType = "2d"
    # Check for type of error to calculate
    elif val[0].lower() == "read":
      if "y" in val[1].lower():
        readFlag = "yes"
  # Last none-check
  if peakNum is None:
    missStr += "\n  Missing peak number."
    errBool = True
  if dataFolders is None:
    missStr += "\n  Missing data folders."
    errBool = True
  if outPath is None:
    missStr += "\n  Missing output path."
    errBool = True

  return errBool, missStr, peakNum, dataFolders, outPath, fitType, readFlag

#---------------------------#---------------------------#
# Takes a string mapped to idices of the fit file and
#  checks to make sure that they are numbers.
#---------------------------#---------------------------#
def ParseIndexString(IdxStr):
  # Split first by commas
  fitLn = IdxStr.split(",")
  errBool, retMsg = False, ""
  FitIdx = set([])
  # Now loop over, and cast as ints and split to ranges.
  for ln in fitLn:
    if "-" in ln:
      try:
        lo, hi = int(ln.split("-")[0]), int(ln.split("-")[1])
        if lo > hi:
          lo, hi = hi, lo
        nl = [x for x in range(lo, hi + 1)]
        for fn in nl:
          FitIdx.add(str(fn))
      except ValueError:
        retMsg += "\n  ERROR: Bad range of fit indices for (%s)\n" % ln
        errBool = True
    else:
      try:
        FitIdx.add(str(ln))
      except ValueError:
        retMsg += "\n  ERROR: Bad fit index for (%s)\n" % ln
        errBool = True
  return errBool, retMsg, [str(x) for x in sorted([int(y) for y in list(FitIdx)])]