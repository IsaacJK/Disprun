#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

# Gen imports
import os, sys, fnmatch, stat, csv
from multiprocessing import Pool
from distutils import dir_util
import subprocess
import time, datetime
# Matpliotlib imports
from matplotlib.backends.backend_pdf import PdfPages
# numpy imports
from numpy import array, asarray, hstack, log, mean, pi, savetxt
# disprun specific libraries
import dr3_gensl as gensl
import dr3_process as prc
import dr3_errors as err
import dr3_plots as dr3plot
import dr3_mathfuncs as mf
# Misc imports
import nmrglue as ng
import pandas as pd
curDir = os.getcwd()
argc = len(sys.argv)

def makeFolder(pathToFolder):
  if not os.path.exists(pathToFolder): 
    os.makedirs(pathToFolder) 

def help():
  print "Disprun v3.06"
  print "Usage is as follows:"
  print " >dr3.py -genpar [output path for input text file]"
  print " >dr3.py -com [input text file]"
  print " >dr3.py -[fit/fit0/fitsl] [input text file]"
  print " >dr3.py -gensl [ouput name for .csv file]"
  print " >dr3.py -clean [folders to be cleaned]"
  print " >dr3.py -swe [R1rho csv] [Error csv] [recombined csv name]"
  print " >dr3.py -cleanvd [input text file]"

  print '''
Disprun v3.06
Usage is as follows:
 >dr3.py -genpar [output path for input text file]
 >dr3.py -com [input text file]
 >dr3.py -[fit/fit0/fitsl] [input text file]
  - fit : fits decaying intensities to monoexponential and props. err
  - fitsl : takes in SLcalibration experiment, fits ints to dec sine func
            outputs fitted cnst12 value in extracted parameters file.
  - fit0 : Takes last delay point and calculates R1p assuming perfect monoexp
 >dr3.py -gensl [ouput name for .csv file]
 >dr3.py -clean [folders to be cleaned]
 >dr3.py -swe [R1rho csv] [Error csv] [recombined csv name]
 >dr3.py -cleanvd [input text file]
'''

def removeFile(filePath):
  subprocess.call(["rm", filePath])

# Checks to see if a given value can be
#  cast as an int.
# Returns true or false
def IsInt(val):
  try:
    int(val)
    return True
  except ValueError:
    return False

# Change directory, run script, return directory
def runCSHScript(parentPath, filePath):
  os.chdir(parentPath)
  subprocess.call(["tcsh", filePath])
  os.chdir(curDir)

def findFiles(directory, keyword):
  foundFiles, foundDirs = [], []
  for root, dirs, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, keyword):
      foundFiles.append(filename)
      foundDirs.append(os.path.join(root, filename))
  return(foundDirs)

if argc == 3 and sys.argv[1].lower() == "-gensl":
  filename = sys.argv[2]
  slpath = os.path.join(curDir, filename)
  gensl.menu(filename, slpath)

####################################
#### Generate an output text file and
#### put in to specified directory
####  Attempts to read folders and assign those automatically
####################################
elif sys.argv[1].lower() == "-genpar" and os.path.isdir(os.path.join(curDir, sys.argv[2])):
  # Path with folders to be renumbered
  dirpath = os.path.join(curDir, sys.argv[2])
  outPath = os.path.join(dirpath, "input_for_dr3.txt")

  ## Try to grab folder ranges to use in the input text file
  # Grab non-numeric folders > 800
  paths = [x for x in os.listdir(dirpath) if IsInt(x) == True and int(x) >= 300]
  ranges = []
  # Sort
  paths = sorted([int(x) for x in paths])

  start, end = 0, 0
  rangeval = "%s-%s"
  for i in range(len(paths)):
    pval,cval,nval = None,None,None
    if i == 0:
      cval = paths[i]
      if len(paths) == 1:
        nval = cval
      else:
        nval = paths[i+1]
    elif i > 0 and i != len(paths) - 1:
      pval = paths[i-1]
      cval = paths[i]
      nval = paths[i+1]
    elif i == len(paths) - 1:
      pval = paths[i-1]
      cval = paths[i]

    if pval is None:
      start = cval
      if cval + 1 != nval:
        end = cval
    elif nval is None:
      end = cval
      if pval + 1 != cval:
        start = cval
    else:
      if pval + 1 != cval and cval + 1 != nval:
        start, end = cval, cval
      elif pval + 1 != cval and cval + 1 == nval:
        start = cval
      elif pval + 1 == cval and cval + 1 != nval:
        end = cval
    
    if start != 0 and end != 0:
      ranges.append(rangeval % (start, end))
      start, end = 0, 0

  if len(ranges) >= 1:
    rangestr = ",".join(ranges)
  else:
    rangestr = "920-1039"

  outInp = '''
# Choose peak number from your peak fitting
PeakNum 1

# '1d' or '2d' for 1D (normal) or pseudo-2D
FitType 2d

# Range of data folders to process
#  Lowest number folder must contain processed fid.com file
DataFolders %s

# Directory to output fits
OutputFolder R1rho_2D

# Do not reprocess data, just re-read fits and re-do error estimation.
Read No

'''
  FILE = open(outPath, "wb")
  FILE.writelines(outInp % rangestr)
  FILE.close()


####################################
#### Fit Peak Ints and Decay Curves ####
####################################
elif (sys.argv[1].lower() == "-fit" or sys.argv[1].lower() == "-fitsl"
      or sys.argv[1].lower() == "-fit0"
      and os.path.isfile(os.path.join(curDir, sys.argv[2]))):
  # Error handling
  errBool, missStr = False, ""
  fit_arg = sys.argv[1].lower()
  # Path to input text file
  inpPath = os.path.join(curDir, sys.argv[2])

  # slpoffs: List to store experimental parameters like SLP, offset, etc
  # r1p: List to store: corr offset, slp (Hz), R1rho, R1rho_err
  slpoffs, dlyints = {}, {}
  r1p_std, r1p_mc, r1p_bs = {}, {}, {}
  r1p_mathematica = {}

  ## Parse input text ##
  # get:
  # -error bool, error string
  # -peak number for peak fitting (i.e. index in test.tab)
  # -sorted (numerical) list of data folders
  # -output path for data
  # -fit type for peak fitting, 1d or 2d
  # -errType is std, mc or bootstrap
  errBool, tstr, peakNum, dataFolders, outPath, fitType, readFlag = prc.ParseInp(inpPath)
  # Check for abortive errors
  missStr += tstr
  err.HandleErrors(errBool, missStr)

  # Make fromatted timestamp
  timeSt = datetime.datetime.fromtimestamp(time.time()).strftime("%Y%m%d%H%M%S")
  # Make output folder and output paths
  outPath = os.path.join(curDir, outPath, timeSt)
  outCopies = os.path.join(outPath, "Copies")
  makeFolder(outPath)
  makeFolder(outCopies)
  # Output folder for extracted expt params
  outParams = os.path.join(outPath, "ExtractedParameters.csv")
  # Output folder for extracted intensities, delays and noise vals
  outDlyInts = os.path.join(outPath, "IntDelayNoise.csv")
  # Output folder for fitted R1rho values
  outR1p_std = os.path.join(outPath, "FittedR1p-std.csv")
  outR1p_mc = os.path.join(outPath, "FittedR1p-mc.csv")
  outR1p_bs = os.path.join(outPath, "FittedR1p-bs.csv")
  outR1p_mathematica = os.path.join(outPath, "FittedR1p-mathematica.csv")
  # Output Mag decay graphs
  outGraphDecMC = os.path.join(outPath, "ExponentialDecays-MC.pdf")
  outGraphDecBS = os.path.join(outPath, "ExponentialDecays-BS.pdf")
  # Output normal MC distriubtion of R1rho
  outMCDist = os.path.join(outPath, "R1pDistribution-MC.pdf")
  # Output Bootstrap distriubtion of R1rho
  outBSDist = os.path.join(outPath, "R1pDistribution-BS.pdf")
  # Copy of input file
  outCPInp = os.path.join(outCopies, sys.argv[2].replace(".txt","")+"-copy.txt")

  # Check that folders and needed files in folder1  exist
  firstPath = os.path.join(curDir, dataFolders[0])
  fld0 = dataFolders[0]

  missStr += tstr
  # Check folders exist, and needed files in folder 0
  errBool, tstr = err.CheckFolders(curDir, dataFolders, fitType, CheckAll=True)
  # Check for abortive errors
  missStr += tstr
  err.HandleErrors(errBool, missStr)

  # If not marked as read-only, then copy all scripts from base directory
  if readFlag == "no":
    # Start copying scripts
    for fn in dataFolders:
      copynames = ["fid.com", "SelR1rho.com", "fit.com", "test.tab"]
      cpTo = os.path.join(curDir, fn)
      prc.CleanDly(cpTo)
      makeFolder(os.path.join(cpTo, "ft"))

      for n in copynames:
        if fn != fld0:
          cpName = os.path.join(firstPath, n)
          subprocess.call(["cp", cpName, cpTo])
        else:
          # If primary directory, make copies of files to copy folder
          cpName = os.path.join(firstPath, n)
          subprocess.call(["cp", cpName, outCopies])

  # Setup PDF Pages for graphing
  plotMCDec = PdfPages(outGraphDecMC)
  plotBSDec = PdfPages(outGraphDecBS)
  plotMCDist = PdfPages(outMCDist)
  plotBSDist = PdfPages(outBSDist)

  # 1. Run scripts
  # 2. Extract experimental parameters
  # 3. Extract delays, noise
  # 4. Fit R1rho, R1rho err.
  for fn in dataFolders:
    runnames = ["fid.com", "SelR1rho.com", "fit.com"]
    runParent = os.path.join(curDir, fn)
    print "\n\n      ~~Processing data for ( %s )~~" % fn
    # If marked as non-read only, then fit peaks
    if readFlag == "no":
      for n in runnames:
        trun = os.path.join(runParent, n)
        runCSHScript(runParent, trun)

    # Grab experimental parameters
    slpoffs[fn] = prc.GrabSLPOffs(runParent, fn)
    
    # Grab new delays, relative intensities, relative noise
    rDIN, nDIN = prc.GrabDlyIntsNoise(runParent, peakNum, fitType)

    # MC-fit monoexp decay
    MCnum = 500 # MC iteration number
    print "\n    -- Fitting monoexponential decay with %s Monte-Carlo & Bootstrap iterations --" % int(MCnum)

    # dly, ints, noiserel = nDIN[:,0], nDIN[:,1], nDIN[:,2]
    dly, ints, noiserel = rDIN[:,0], rDIN[:,1], rDIN[:,2]
    if fit_arg == "-fit0":
      dly = dly[:-1]
      ints = ints[:-1]
      noiserel = noiserel[:-1]
    # Put delays and ints to a dict for later export
    dlyints[fn] = array(array([rDIN[:,0], rDIN[:,1], rDIN[:,2],
                                     nDIN[:,1], nDIN[:,2]]).T)


    # Grab some mag and delay values
    Tmax = dly.max()
    magMin = ints.min()
    magMax = ints.max()
    
    # If peak intensity drops below zero, then rescale to positive
    if magMin <= 0.:
      print " !!! WARNING intensity of peak goes below 0 !!!"
      magMin = 0.1 - magMin
      magMax = magMax + magMin

    # R1rho fitted values and errors and distributions
    R1p_mu_std, R1p_mu_MC, R1p_mu_BS = None, None, None
    R1p_sigma_std, R1p_sigma_MC, R1p_sigma_BS = None, None, None
    # N-number of Prexponential factors and R1rhos
    FitVals_MC, FitVals_BS = None, None

    if fit_arg == "-fit" or fit_arg == "-fit0":
      # Initial guess for fit
      R1p_t = -1./Tmax*log(magMin/magMax)
      p0 = array([max(ints), R1p_t])
      func = mf.expdec

    elif sys.argv[1].lower() == "-fitsl":
      # Initial guess for damped sine fit:
      #  - Amplitude
      #  - R1
      #  - Omega (~cnst12)
      #  - Phi (small, unknown)
      omg = float(slpoffs[fn]['cnst12']) * 2. * pi
      p0 = array([max(ints), 10., omg, 1.])
      func = mf.dampedsine

    # Calculate R1rho mean and sigma from MC distribution
    Std_Fit, Std_Err = mf.StandardError(func, p0, ints, dly, noiserel)
    A_std, R1p_mu_std = Std_Fit[0], Std_Fit[1]
    R1p_sigma_std = Std_Err[1]
    # # Calculate R1rho mean and sigma from MC distribution
    MC_Fit, MC_Err, FitVals_MC = mf.MonteCarlo(func, p0, ints,
                                                dly, noiserel, MCnum, StdErr = False)
    R1p_mu_MC, R1p_sigma_MC = MC_Fit[1], MC_Err[1]
    
    if sys.argv[1].lower() == "-fit":
      # Calculate R1rho mean and sigma from Bootstrap approach
      BS_Fit, BS_Err, FitVals_BS = mf.Bootstrap(func, p0, ints,
                                                  dly, noiserel, MCnum, StdErr = False)
      R1p_mu_BS, R1p_sigma_BS = BS_Fit[1], BS_Err[1]

    elif sys.argv[1].lower() == "-fitsl":
      slpoffs[fn]['cnst12-cal_MC'] = MC_Fit[2] / (2. * pi)
      slpoffs[fn]['cnst12-calErr_MC'] = MC_Err[2] / (2. * pi)

    # Calculate reduce-chi squares for the different fits
    redchi2_std = mf.redchi2(func, Std_Fit, dly, ints, noiserel)
    redchi2_MC = mf.redchi2(func, MC_Fit, dly, ints, noiserel)
    if sys.argv[1].lower() == "-fit":
      redchi2_BS = mf.redchi2(func, BS_Fit, dly, ints, noiserel)
    # Generate r1p dictionary with fitted R1rho values and errors
    r1p_std[fn] = [fn, slpoffs[fn]['cnst28i'], slpoffs[fn]['cnst12'],
                   R1p_mu_std, R1p_sigma_std, redchi2_std]

    # Monte-Carlo Error
    # - R1rho = weighted best-fit
    # - R1rho_error = MC Error
    #   Since R1rho from weighted best-fit, red. chi^2 is from std
    r1p_mc[fn] = [fn, slpoffs[fn]['cnst28i'], slpoffs[fn]['cnst12'],
                   R1p_mu_std, R1p_sigma_MC, redchi2_std]

    if sys.argv[1].lower() == "-fit":
      # Bootstrap Error
      # - R1rho = weighted best-fit
      # - R1rho_error = Bootstrap Error
      #   Since R1rho from weighted best-fit, red. chi^2 is from std
      r1p_bs[fn] = [fn, slpoffs[fn]['cnst28i'], slpoffs[fn]['cnst12'],
                     R1p_mu_std, R1p_sigma_BS, redchi2_std]

    # Similar to mathematica fits: MC mean R1rho, Std error
    r1p_mathematica[fn] = [fn, slpoffs[fn]['cnst28i'], slpoffs[fn]['cnst12'],
                   R1p_mu_MC, R1p_sigma_std, redchi2_MC]

    # Plot MC distribution of fitted R1rhos
    plotMCDist.savefig(dr3plot.PlotDist(FitVals_MC, R1p_mu_MC, R1p_sigma_MC,
                                        slpoffs[fn]['cnst12'], slpoffs[fn]['cnst28i'],
                                        slpoffs[fn]['FN']), index=1)

    if sys.argv[1].lower() == "-fit":
      # Plot BS distribution of fitted R1rhos
      plotBSDist.savefig(dr3plot.PlotDist(FitVals_BS, R1p_mu_BS, R1p_sigma_BS,
                                          slpoffs[fn]['cnst12'], slpoffs[fn]['cnst28i'],
                                          slpoffs[fn]['FN']), index=1)

    # Plot monoexponential decay curve
    plotMCDec.savefig(dr3plot.PlotDecay(dly, ints, noiserel, FitVals_MC, R1p_mu_MC, R1p_sigma_MC,
                                      slpoffs[fn]['cnst12'], slpoffs[fn]['cnst28i'],
                                      slpoffs[fn]['FN'], redchi2_MC, MC_Fit, func))
    if sys.argv[1].lower() == "-fit":
      plotBSDec.savefig(dr3plot.PlotDecay(dly, ints, noiserel, FitVals_BS, R1p_mu_BS, R1p_sigma_BS,
                                        slpoffs[fn]['cnst12'], slpoffs[fn]['cnst28i'],
                                        slpoffs[fn]['FN'], redchi2_BS, BS_Fit, func))
    print " -----------------------------"
  
  # Close plotting function
  plotMCDec.close()
  plotBSDec.close()
  plotMCDist.close()
  plotBSDist.close()
  # Write out the extracted experimental parameters
  FILE = open(outParams, "wb")
  FILE.write(",".join([str(x) for x in slpoffs[dataFolders[0]].keys()]) + "\n")
  for fn in dataFolders:
    FILE.write(",".join([str(x) for x in slpoffs[fn].values()]) + "\n")
  FILE.close()

  # Write out the MC fitted R1rho values and error
  FILE = open(outR1p_std, "wb")
  FILE2 = open(outR1p_mc, "wb")
  FILE3 = open(outR1p_bs, "wb")
  FILE4 = open(outR1p_mathematica, "wb")
  FILE.write("Folder, Offset, SLP, R1p_STD, R1p_err_STD, red_chi_sq\n")
  FILE2.write("Folder, Offset, SLP, R1p_MC, R1p_err_MC, red_chi_sq\n")
  FILE3.write("Folder, Offset, SLP, R1p_BS, R1p_err_BS, red_chi_sq\n")
  FILE4.write("Folder, Offset, SLP, R1p_MC, R1p_err_STD, red_chi_sq\n")
  for fn in dataFolders:
    FILE.write(",".join([str(x) for x in r1p_std[fn]]) + "\n")
    FILE2.write(",".join([str(x) for x in r1p_mc[fn]]) + "\n")
    if sys.argv[1].lower() == "-fit": 
      FILE3.write(",".join([str(x) for x in r1p_bs[fn]]) + "\n")  
    FILE4.write(",".join([str(x) for x in r1p_mathematica[fn]]) + "\n")  
  FILE.close()
  FILE2.close()
  FILE3.close()
  FILE4.close()

  # # Write out the delays and intensities
  # FILE = open(outDlyInts, "wb")
  # FILE.write("Folder Number, Delay (sec), Raw Intensity, Raw Noise, Normalized Intensity, Normalized Noise\n")
  # for fn in dataFolders:
  #   for din in [x for x in dlyints[fn]]:
  #     FILE.write(str(fn) + ",")
  #     FILE.write(",".join([str(y) for y in din]) + "\n")
  # FILE.close()
  mdf = pd.DataFrame()
  for fn in dataFolders:
    # Get numpy array of corrected offset and slp to match number of dlys
    tsloff = array([[fn, slpoffs[fn]['cnst28i'], slpoffs[fn]['cnst12']]
                     for _ in range(dlyints[fn].shape[0])])
    # Stack array of offsets/slps to array of dlys/intensities
    intsout = hstack((tsloff, dlyints[fn]))
    # Column headers
    intsout_h = ["Folder", "Offset", "SLP", "Dly",
                 "Int", "Int_err", "nInt", "nInt_err"]
    # append dataframe of this array to master dataframe to be written out
    mdf = mdf.append(pd.DataFrame(intsout, columns=intsout_h))
  # reset index of dataframe
  mdf = mdf.reset_index(drop=True)
  mdf.to_csv(outDlyInts, sep=",", index=False)
  # Make copies of fit files for archive purposes
  subprocess.call(["cp", inpPath, outCPInp])

####################################
#### Grab cnst30 and o3p        ####
####################################
elif (argc == 3 and sys.argv[1].lower() == "-cs"):
  pathC = os.path.join(curDir, sys.argv[2], "1420")
  pathN = os.path.join(curDir, sys.argv[2], "920")
  acqusC = os.path.join(curDir, sys.argv[2], "1420", "acqus")
  acqusN = os.path.join(curDir, sys.argv[2], "920", "acqus")
  C_H, N_H, o2p, o3p = None, None, None, None
  if os.path.isfile(acqusC):
    vals = prc.GrabSLPOffs(pathC, "1420")
  if os.path.isfile(acqusN):
    vals = prc.GrabSLPOffs(pathN, "920")
    print vals['cnst30']
    print vals['o3p']

####################################
#### Swap errors                ####
#### dr3.py -swaperr r1p.csv r1perr.csv newfile.csv
####################################
elif (4 <= argc <= 5 and sys.argv[1].lower() == "-swe" 
      and os.path.isfile(os.path.join(curDir, sys.argv[2]))
      and os.path.isfile(os.path.join(curDir, sys.argv[3]))):
  # output directory for recombined error
  if argc == 5:
    if sys.argv[4].endswith(".csv"):
      outpath = os.path.join(curDir, sys.argv[4])
    else:
      outpath = os.path.join(curDir, sys.argv[4] + ".csv")
  else:
    outpath = os.path.join(curDir, "SwappedError.csv", "wb")

  # Path for R1rho values and R1rho error
  r1ppath = os.path.join(curDir, sys.argv[2])
  errpath = os.path.join(curDir, sys.argv[3])

  # Open and extract data from input files
  FILE = open(r1ppath, "rU")
  r1pd = [x.strip().split(",") for x in FILE]
  FILE.close()

  # Make array
  r1pd = array(r1pd)

  # Open and extract data from input files
  FILE = open(errpath, "rU")
  errd = [x.strip().split(",") for x in FILE]
  FILE.close()

  # Make array
  errd = array(errd)

  # Swap error columns
  r1pd[:,4] = errd[:,4]
  
  FILE = open(outpath, "wb")
  for line in r1pd:
    FILE.write(",".join(line) + "\n")
  FILE.close()
####################################
#### Replace vdlist 'ms'        ####
####################################
elif (3 <= argc <= 4 and sys.argv[1].lower() == "-cleanvd" 
      and os.path.isfile(os.path.join(curDir, sys.argv[2]))):
  # Error handling
  errBool, missStr = False, ""

  # Path to input text file
  inpPath = os.path.join(curDir, sys.argv[2])

  ## Parse input text ##
  # get:
  # -error bool, error string
  # -peak number for peak fitting (i.e. index in test.tab)
  # -sorted (numerical) list of data folders
  # -output path for data
  # -fit type for peak fitting, 1d or 2d
  # -readFlag : 'yes' or 'no', reads data only, does not re-fit peaks
  errBool, tstr, peakNum, dataFolders, outPath, fitType, readFlag = prc.ParseInp(inpPath)
  # Check for abortive errors
  missStr += tstr
  err.HandleErrors(errBool, missStr)

  # Go through data folders and fix vdlists
  for p in dataFolders:
    tpath = os.path.join(curDir, p)
    prc.CleanDly(tpath)
    # # Get path to vdlist in first directory
    # vdlistpath = os.path.join(tpath, "vdlist")

    # if os.path.isfile(vdlistpath):
    #   FILE = open(vdlistpath, "rU")
    #   delays = [x.strip() for x in FILE]
    #   FILE.close()

    #   if "ms" in delays[0]:
    #     delays = [str(float(x.replace("ms", ""))/1e3) for x in delays if len(x.replace("ms", "")) > 0]
    #     FILE = open(vdlistpath, "wb")
    #     for line in delays:
    #       FILE.write(line + "\n")
    #     FILE.close()
    # else:
    #   print "Path for %s does not exist." % p
  
####################################
#### Generate SelR1rho.com file ####
####################################
elif (3 <= argc <= 4 and sys.argv[1].lower() == "-com" 
      and os.path.isfile(os.path.join(curDir, sys.argv[2]))):
  # Error handling
  errBool, missStr = False, ""

  # Path to input text file
  inpPath = os.path.join(curDir, sys.argv[2])

  ## Parse input text ##
  # get:
  # -error bool, error string
  # -peak number for peak fitting (i.e. index in test.tab)
  # -sorted (numerical) list of data folders
  # -output path for data
  # -fit type for peak fitting, 1d or 2d
  # -readFlag : 'yes' or 'no', reads data only, does not re-fit peaks
  errBool, tstr, peakNum, dataFolders, outPath, fitType, readFlag = prc.ParseInp(inpPath)
  # Check for abortive errors
  missStr += tstr
  err.HandleErrors(errBool, missStr)

  # Grab first path in list
  firstPath = os.path.join(curDir, dataFolders[0])
  # Clean vdlist in first path, just to be sure
  prc.CleanDly(firstPath)
  # Get path to vdlist in first directory
  vdlistpath = os.path.join(firstPath, "vdlist")

  # Set phase correction, if given
  p0 = 0.0
  if argc == 4:
    try:
      p0 = float(sys.argv[3])
    except ValueError:
      print "Invalid p0 given, setting as '0.0'"
      p0 = 0.0

  # Check vdlist exists
  if not os.path.isfile(vdlistpath):
    print "No vdlist found."
    sys.exit(-1)

  # Read in delays
  FILE = open(vdlistpath, "rU")
  delays = [x.strip() for x in FILE]
  FILE.close()

  # Check to make sure ms not in vdlist
  if "ms" in delays[0]:
    delays = [float(x.replace("ms", ""))/1e3 for x in delays]

  # Create SelR1rho.com file in path specified
  selPath = prc.WriteSelCom(firstPath, delays, fitType, p0)

  # Change permissions of com script
  # read,write,execute by owner
  os.chmod(selPath, stat.S_IRWXU)

  # Check to see if the output of the fid.com file
  #  already exists so that we can run SelR1rho.com
  if os.path.isfile(os.path.join(firstPath, "test.fid")):
    # Run SelR1rho.com
    runCSHScript(firstPath, selPath)

###################################
#### Clean Folders Recursively ####
###################################
elif argc == 3 and sys.argv[1].lower() == "-clean" \
     and os.path.isdir(os.path.join(curDir, sys.argv[2])):
  delfiles = ["R1rho1D*", "Fit-*","RUN","DOIT","COPIER","COMBINE","out.ft2","sim.ft1",
              "test.fid", "copySel.com", "BrukParams.txt","autoFit.com","axt.tab",
              "dif.ft1","pk.tcl","pkROI_1d.dat", "test*.ft*", "dif.ft1", "test.tab",
              "tau.list", "Fit-*.tab", "noise-*.txt", "R1rho1D-*.ft2", "SelR1r-ints_*.txt",
              "selectedT1fits.txt", "Extracted_Parameters.csv", "SelR1*", "fid.com",
              "fit.tab", "fit.com"]
  delMeDir = []
  if sys.argv[2] != ".":
    cleanFolder = os.path.join(curDir, sys.argv[2])
  else:
    cleanFolder = curDir
  totalSize = 0.0
  pool = Pool(processes=8)

  for i in delfiles:
    delPaths = findFiles(cleanFolder, i)
    for x in delPaths:
      print x
      tpathnum = None
      try:
        tpathnum = int(x.split("/")[-2])
      except ValueError:
        try:
          tpathnum = int(x.split("/")[-3])
        except ValueError:
          if x.split("/")[-2] == "OnRes" or x.split("/")[-3] == "OnRes":
            tpathnum == "OnRes"
          else:
            print "Not valid path ( %s )" % x.split("/")[-2]
      if 100 <= tpathnum <= 3000 or tpathnum == "OnRes":
        totalSize += os.stat(x).st_size
        # delMeDir.append(x)
        subprocess.call(["rm",x])
        # removeFile(x)
  # results = [pool.apply_async(removeFile,args=(x,)) for x in delMeDir]
  # output = [p.get() for p in results]
  print "Removed %s Mb of files from '%s'" % (round(totalSize/1e6,1), cleanFolder)
  
###################################
#### Estimate S/N from fit.tab ####
###################################
elif (3 <= argc <= 4 and sys.argv[1].lower() == "-sn"
     and os.path.isfile(os.path.join(curDir, sys.argv[2]))):
  
  # Path to tab file with signal and noise
  tabpath = os.path.join(curDir, sys.argv[2])

  # NMRglue read tab
  pc,pf,rec = ng.fileio.pipe.read_table(tabpath)

  # Signal dictionary
  sigdict = {}

  # Find spectral noise
  for i in pc:
    tline = i.replace(",", "").split(" ")
    if "Noise" in tline[1]:
      noise = float(tline[2])

  print "\nNoise is: %.0f" % noise
  # Get signal values
  for i in rec:
    print "  Peak #%s: S/N ( %.2f ), Signal ( %.0f )" % (i[0], float(i[9])/noise, i[9])
  print ""

###################################
#### Generate R1rho folders and setup scripts
###################################
elif  (sys.argv[1].upper() == "-15N" or sys.argv[1].upper() == "-13C"
    and os.path.isfile(os.path.join(curDir, sys.argv[2]))
    and os.path.isdir(os.path.join(curDir, sys.argv[3]))):
  data = []
  nuclei = sys.argv[1].upper()
  csvpath = os.path.join(curDir, sys.argv[2])
  folderpath = os.path.join(curDir, sys.argv[3])
  if argc == 5:
    fldname = sys.argv[4]
  else:
    if "15N" in sys.argv[1].upper():
      fldname = "15N"
    elif "13C" in sys.argv[1].upper():
      fldname = "13C"
  if sys.argv[3][-1:] == "/":
    foldername = sys.argv[3][:-1]
  else:
    foldername = sys.argv[3]

  if os.path.isfile(csvpath) == True and os.path.exists(folderpath) == True:
    with open(csvpath, "rU") as csvfile:
      readme = csv.reader(csvfile, delimiter=",")
      for row in readme:
        data.append(row)
    # Check for a header in the file
    try:
      float(data[0][0])
    # If can't cast a float, assume header and delete first row
    except ValueError:
      del(data[0])

    if nuclei == "-13C":
      setparName = "%s-OffResSetpar" % fldname
    elif nuclei == "-15N":
      setparName = "%s-OffResSetpar" % fldname

    FILE = open(setparName, "wb")
    FILE.write("re %s %s 1\n\n" % (fldname, data[0][0]))
    for i in data:
        #Folder #
      FILE.write("re " + str(i[0]) + " " + "1\n")
       #UNCORRECTED Offset frequency, Hz
      FILE.write("cnst28 " + str(i[1]) + "\n")
      FILE.write("%s-setpar\n" % fldname)
        #Spinlock power, Hz
      FILE.write("cnst12 " + str(i[3]) + "\n")
        #Spinlock power, dB
      FILE.write("pldb23 " + str(i[4]) + "\n")
      if nuclei == "-13C":
        if float(i[1]) != 0.:
          FILE.write("zgoptns -DC13_OFFSET\n")
        else:
          FILE.write("\n")
        
      elif nuclei == "-15N":
        if float(i[1]) != 0.:
          FILE.write("zgoptns -DN15_OFFSET\n")
        else:
          FILE.write("\n")
      FILE.write("\n")
    FILE.close()

    if nuclei == "-13C":
      FILE = open("%s-run" % fldname, "wb")
    elif nuclei == "-15N":
      FILE = open("%s-run" % fldname, "wb")
    FILE.write("re %s %s 1\n\n" % (fldname, data[0][0]))
    for i in data:
        #Folder number
      FILE.write("re " + str(i[0]) + " " + "1\n")
      FILE.write("zg\n")
      FILE.write("\n")
    FILE.close()

    for i in data:
      dir_util.copy_tree(folderpath, str(i[0]))

    FILE = open(os.path.join(curDir, "%s-setpar" % fldname), "wb")
    if nuclei == "-13C":
      FILE.write('''
o1p 4.7
p1 11
te 298.2
d1 1.5
ns 16
ds 2
rg 203

cnst30 7
o2p 159
o3p 190

spdb1 38
spdb2 28
p11 4600

cnst11 9000
pldb24 3.2

d31 1ms
d30 2ms
vdlist 13C_zero
nbl 1
1 td 1

pcpd2 105
pldb12 -0.82
pcpd3 220
pldb16 -3.9
''')

    elif nuclei == "-15N":
      FILE.write('''
o1p 4.7
p1 11
te 298.2
d1 1.5
ns 16
ds 2
rg 203

cnst30 7
o2p 158
o3p 145

cnst11 6900
pldb24 0.0
p25 1600
spdb5 29.5

spdb1 38
spdb2 18
p11 10500

d31 1ms
d30 2ms
vdlist 15N_zero
nbl 1
1 td 1

pcpd2 120
pldb12 1000
pcpd3 200
pldb16 -4
''')
    FILE.close()
  else:
    print ".csv file or folder does not exists."

###################################
#### Generate R1rho folders and setup scripts
###################################
elif  (sys.argv[1].lower() == "-ints"
      and os.path.isfile(os.path.join(curDir, sys.argv[2]))
      and os.path.isfile(os.path.join(curDir, sys.argv[3]))):
    pars_dir = os.path.join(curDir, sys.argv[2])
    ints_dir = os.path.join(curDir, sys.argv[3])
    p_df = pd.read_csv(pars_dir, sep=",")
    i_df = pd.read_csv(ints_dir, sep=",")
    i_df = i_df.rename(columns = {"Folder Number": "FN",
                                  " Delays (sec)": "dly",
                                  " Raw Intensity": "ri",
                                  " Raw Noise": "rn",
                                  " Normalized Intensity": "ni",
                                  " Normalized Noise": "nn"})
    o_df = pd.merge(i_df, p_df, on="FN")
    o_df.to_csv("test.csv", sep=",")
    # print p_df['FN']
    # print i_df['FN']
    # print i_df.columns
else: help()

