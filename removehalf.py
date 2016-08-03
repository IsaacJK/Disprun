import os, sys, fnmatch

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
  print " >removehalf.py [csv file] [+/-]"


def findFiles(directory, keyword):
  foundFiles, foundDirs = [], []
  for root, dirs, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, keyword):
      foundFiles.append(filename)
      foundDirs.append(os.path.join(root, filename))
  return(foundFiles, foundDirs)

if argc == 3 and fileExists(os.path.join(curDir, sys.argv[1])):
  indata,outdata = [],[]
  inpath = os.path.join(curDir, sys.argv[1])
  sign = sys.argv[2]
  FILE=open(inpath, "rU")
  indata = [x.strip().split(",") for x in FILE]
  FILE.close
  ct = 1
  for idx in range(len(indata)):
    toff = int(indata[idx][1])    
    if sign == "+":
      if toff > 0:
        ct += 1
        if ct % 2 == 0:
         None
        else:
          outdata.append(indata[idx])
      else:
        outdata.append(indata[idx])
    elif sign == "-":
      if toff < 0:
        ct += 1
        if ct % 2 == 0:
         None
        else:
          outdata.append(indata[idx])
      else:
        outdata.append(indata[idx])

  for idx,line in enumerate(outdata):
    line[0] = str(idx+920)

  FILE = open("outTrimFold.csv", "wb")
  for line in outdata:
    FILE.write(",".join(line))
    FILE.write("\n")
  FILE.close

    # if int(line[1]) >

else: help()
