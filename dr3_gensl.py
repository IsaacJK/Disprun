#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

import os, sys, fnmatch
import dr3_SLEqns as SLEqns

curDir = os.getcwd()
argc = len(sys.argv)

def fileExists(pathToFile):
  try: 
    with open(pathToFile) as f: return(True)
  except IOError as e: return(False) 

def makeFolder(pathToFolder):
  if not os.path.exists(pathToFolder): 
    os.makedirs(pathToFolder) 

def help():
  print "Usage is as follows:"
  print " >something.py [otherstuff]"

def IsType(s,t):
    if t == "Float":
      try: 
          float(s)
          return True
      except ValueError:
          return False
    elif t == "Int":
      try: 
          int(s)
          return True
      except ValueError:
          return False

def findFiles(directory, keyword):
  foundFiles, foundDirs = [], []
  for root, dirs, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, keyword):
      foundFiles.append(filename)
      foundDirs.append(os.path.join(root, filename))
  return(foundFiles, foundDirs)

def menu(filename, fullpath):
  MenuLoop, MenuLoop2 = True, True
  AltLoop = True
  CalSet = False
  OnRes = []
  OffRes = []
  OutData = []
  nuclei, bf, salt = None,None,None
  # For offresonance
  spinlocks = []
  vals = 24
  multiplier = 3.5    

  while MenuLoop == True:
    AltLoop = True
    MenuLoop2 = True
    print "Generating %s R1p points (%s OnRes & %s OffRes)" % ((len(OnRes)+len(OffRes)), len(OnRes), len(OffRes))
    print "Using %s OffRes SLPs and %s offsets each SLP." % (len(spinlocks),vals)
    print " Menu Options:"
    print " (1) [OnRes] Spinlock Powers/Points"
    print " (2) [OffRes] Spinlock Powers/Offsets"
    print " (3) Calibration curve"
    print " (4) Reset"
    print " (5) Write out list"
    print " (0) Quit"

    mOpt = raw_input(">")
    print ""

    ## Quit ##
    if mOpt == "0":
      MenuLoop = False

    ## Set On-Res Points ##
    elif mOpt == "1":
      onLoop = True
      while AltLoop == True:
        print "Choose On-Resonance Profiles"
        print " (1) Manually add On-Res Points"
        print " (2) Add Standard 13C On-Res Points"
        print " (3) Add Standard 15N On-Res Points"
        print " (4) Remove On-Res Points"
        print " (5) Save and return"
        
        aOpt = raw_input(">")
        print ""

        # Manually add on-resonance data points
        if aOpt == "1":
          # Temp on-res SL list
          tOR = []
          while onLoop == True:
            if len(tOR) != 0:
              print "Spinlock Powers: ", tOR
            print "Keep adding spinlock powers. Type 'quit' to stop."
            print "Add '-' before spinlock power to remove it."
            tSL = raw_input(">")
            
            if IsType(tSL, "Float") == True and tSL[0] != "-":
              if float(tSL) > 3500.:
                print "Value not added, spinlock power (%s) is too high." % tSL
              else:
                tOR.append(tSL)

            elif tSL[0] == "-":
              tSL = tSL[1:]
              tOR = [x for x in tOR if x != tSL]

            elif tSL.lower() == "quit":
              for i in tOR:
                OnRes.append([0.,0.,i])
              onLoop = False

        elif aOpt == "2":
          tSL = [150.,200.,250.,300.,350.,400.,500.,600.,700.,800.,900.,
                 1000.,1200.,1400.,1600.,1800.,2000.,2500.,3000.,3500.]
          for i in tSL:
            OnRes.append([0.,0.,i])
          AltLoop = False

        elif aOpt == "3":
          tSL = [100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,
                 700.,800.,900.,1000.,1200.,1400.,1600.,1800.,2000.]
          for i in tSL:
            OnRes.append([0.,0.,i])
          AltLoop = False

        # Reset on-res points
        elif aOpt == "4":
          OnRes = []

        # Save and return
        elif aOpt == "5":
          AltLoop = False
        else:
          print "Invalid option. Try again."

    ## Set Off-Res Points ##
    elif mOpt == "2":
      while MenuLoop2 == True:
        AltLoop = True
        if len(spinlocks) != 0:
          print "Spinlock Powers: ", spinlocks
          print "Number OffRes R1p points (%s)" % (vals * len(spinlocks))
        else:
          print "No spinlock powers given for off-resonance data."
        print "Offset Multiplier (%s) - Number Offsets per SLP (%s)" % (multiplier,vals)
        print ""

        print " Offset Options:"
        print " (1) Offset Multiplier"
        print " (2) Spinlock offsets per spinlock power"
        print " (3) Add/Remove spinlock powers"
        print " (4) Reset"
        print " (5) Save and return"
        print " (0) Return without saving"

        aOpt = raw_input(">")
        print ""

        if aOpt == "0":
          MenuLoop2 = False

        elif aOpt == "1":
          while AltLoop == True:
            print "Define offset multiplier (numerical value)"
            multiplier = raw_input(">")

            if IsType(multiplier, "Float") == True:
              multiplier = float(multiplier)
              AltLoop = False
            else:
              print "Define a numerical value\n"

        elif aOpt == "2":
          while AltLoop == True:
            print "Define number of offsets per spinlock power"
            vals = raw_input(">")

            if IsType(vals, "Int") == True:
              vals = int(vals)
              AltLoop = False
            else:
              print "Define a numerical value\n"

        elif aOpt == "3":
          while AltLoop == True:
            if len(spinlocks) != 0:
              print "Spinlock Powers: ", spinlocks
            print "Keep adding spinlock powers. Type 'quit' to stop."
            print "Add '-' before spinlock power to remove it."
            tSL = raw_input(">")
            
            if IsType(tSL, "Float") == True and tSL[0] != "-":
              if float(tSL) > 3500.:
                print "Value not added, spinlock power (%s) is too high." % tSL
              elif tSL in spinlocks:
                print "Spinlock is already in list."
              else:
                spinlocks.append(tSL)

            elif tSL[0] == "-":
              tSL = tSL[1:]
              spinlocks = [x for x in spinlocks if x != tSL]

            elif tSL.lower() == "quit":
              AltLoop = False

        elif aOpt == "4":
          spinlocks = []
          OffRes = []

        elif aOpt == "5":
          for i in spinlocks:
            SL = float(i)
            maxoff = multiplier*SL
            newval, minval = 0.,10.
            #Check for non-even vals
            if vals % 2 == 0:
              incr = round(maxoff/((vals/2)-1))
              # Positive offsets
              OffRes.append([minval,-minval,SL])
              for y in range(0,(vals/2)-1):
                newval += incr
                OffRes.append([newval,-newval,SL])
              # Negative offsets
              newval = 0.
              OffRes.append([-minval,minval,SL])
              for y in range(0,(vals/2)-1):
                newval += incr
                OffRes.append([-newval,newval,SL])

            else:
              incr = round(maxoff/(((vals-1)/2)-1))

              # Positive offsets
              OffRes.append([0.,0.,SL])
              OffRes.append([minval,-minval,SL])
              for y in range(0,((vals-1)/2)-1):
                newval += incr
                OffRes.append([newval,-newval,SL])
              # Negative offsets
              newval = 0.
              OffRes.append([-minval,minval,SL])
              for y in range(0,((vals-1)/2)-1):
                newval += incr
                OffRes.append([-newval,newval,SL])
          MenuLoop2 = False

    ## Set Calibration Curve ##
    elif mOpt == "3":
      while AltLoop == True:
        print "Select Nuclei (C or N)"
        nuclei = raw_input(">")
        print ""

        if nuclei.lower() == "c" or nuclei.lower() == "n":
          AltLoop = False
        else:
          print "Input correct nuclei type (C or N)."
      
      AltLoop = True

      while AltLoop == True:
        print "Select Spectrometer frequency (600, 700, or 700ss)"
        bf = raw_input(">")
        print ""

        if bf == "600" or bf == "700" or bf == "700ss":
          AltLoop = False

        else:
          print "Input spectrometer frequency (600, 700, or 700ss)."

      AltLoop = True

      while AltLoop == True:
        print "Select low or high salt (low or high)"
        salt = raw_input(">")
        print ""

        if salt.lower() == "low" or salt.lower() == "lo":
          salt = "low"
          AltLoop = False
        elif salt.lower() == "high" or salt.lower() == "hi":
          salt = "high"
          AltLoop = False
        else:
          print "Input salt concentration (low or high)."          
      CalSet = True

    # Reset on-res and off-res SLPs
    elif mOpt == "4":
      OnRes = []
      OffRes = []

    # Write out values
    elif mOpt == "5":
      if CalSet == True:
        if nuclei == "N":
          start = 900
        else:
          start = 1400

        end = None

        if len(OnRes) != 0:
          for n,i in enumerate(OnRes,start):
            i2 = [str(x) for x in i]
            OutData.append([str(n)] + i2 + [str(SLEqns.calcPLDB23(i[2],bf,nuclei,salt))])
            end = n+1
        else:
          end = start

        for n,i in enumerate(OffRes,end):
          i2 = [str(x) for x in i]
          OutData.append([str(n)] + i2 + [str(SLEqns.calcPLDB23(i[2],bf,nuclei,salt))])
        FILE = open(fullpath, "wb")
        FILE.write("Folder,UncorrOffset,CorrOffset,SLP(Hz),SLP(dB)\n")
        for i in OutData:
          FILE.write(str(",").join(i)+"\n")
        FILE.close()
      else:
        print "  ERROR: Calibration curve not defined!"

  # return((filename,fullpath))
