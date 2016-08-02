#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

import math

  # Returns a PLDB23 float with 2 decimal places, IE 17.60,
  #  dependening on the Hz value given, calculated from the
  #  corresponding calibration curve equation.
def calcPLDB23(HzVal=None, freq=None, nuclei=None, salt=None):
  pldb23 = 1000.0
  if freq == '700ss' and nuclei == 'N' and salt == "low":
    # Construct: wtTAR
    # Temperature: 298K
    # pH: 6.4
    # Salt Conc: 25 mM NaCl
    # Phosphate Conc: 15 mM PO4
    # Resonance: rU42-N3
    # Spectrometer: Bruker 700s
    # Seq type: Sel 15N R1rho
    # Number of points: 30 x 2
    # SL Cal Range: 50-3000 Hz
    # Error in duplicates: No
    # Date collected: 08-08-2016
    # Originator: Isaac K
    pldb23 = (-8.782*math.log(float(HzVal))+52.223)

  elif freq == '700' and nuclei == 'C' and salt == "low":
    # Construct: f.lbl hpTG-GTG
    # Temperature: 298K
    # pH: 5.88
    # Salt Conc: 25 mM NaCl
    # Phosphate Conc: 15 mM PO4
    # Resonance: A3-C8
    # Spectrometer: Bruker 700
    # Seq type: Presat 13C Qi mod
    # Number of points: 30 x 1
    # SL Cal Range: 50-3600 Hz
    # Error in duplicates: No
    # Date collected: 10-15-2015
    # Originator: Isaac K
    pldb23 = (-8.74883*math.log(float(HzVal))+67.62061)

    ## Old Cal Below ##
    # # Construct: A-site
    # # Temperature: 298K
    # # pH: 6.4
    # # Salt Conc: 25 mM NaCl
    # # Phosphate Conc: 15 mM PO4
    # # Resonance: A93 C2
    # # Spectrometer: Bruker 700
    # # Number of points: 20
    # # SL Cal Range: 150-3500 Hz
    # # Error in duplicates: None
    # # Date collected: ???
    # # Originator: Isaac K
    # pldb23 = (-8.866*math.log(float(HzVal))+68.481)

  elif freq == '700' and nuclei == 'C' and salt == "high":
    # Construct: DickersonGT G/T.lbl
    # Temperature: 298K
    # pH: 5.93
    # Salt Conc: 125 mM NaCl
    # Phosphate Conc: 15 mM PO4
    # Resonance: G4-C8 & T7-C6, Avg. stdev
    # Spectrometer: Bruker 700
    # Number of points: 25x2
    # SL Cal Range: 50-3500 Hz (did not use T7-C6 50/100Hz SLP)
    # Error in duplicates: Stdev
    # Date collected: 03-26-2015
    # p1: 13.63
    # Fit type: Origin weighted fit
    # Originator: Isaac K
    pldb23 = (-8.7579*math.log(float(HzVal))+67.55362)

  elif freq == '600' and nuclei == 'C' and salt == "low":
    # Construct: hpCG-GCG f.lbl
    # Temperature: 298K
    # pH: 5.88
    # Salt Conc: 25 mM NaCl
    # Phosphate Conc: 15 mM
    # Resonance: G4-C8 and T17-C6
    # Spectrometer: Bruker 600
    # Number of points: 30x2
    # SL Cal Range: 50-3600 Hz
    # Error in duplicates: None
    # Date collected: 04-11-2015
    # Originator: Isaac K
    pldb23 = (-8.76405*math.log(float(HzVal))+69.35298)

  elif freq == '700' and nuclei == 'N' and salt == "low":         
    # Construct: f.lbl hpTG-GTG
    # Temperature: 298K
    # pH: 6.0
    # Salt Conc: 25 mM NaCl
    # Phosphate Conc: 15 mM
    # Resonance: hpTG T17-N3
    # Spectrometer: Bruker 700
    # Number of points: 25 x1
    # SL Cal Range: 50-2000 Hz
    # Error in duplicates: No
    # Date collected: 10-15-2015
    # Originator: Isaac K
    pldb23 = (-8.74013*math.log(float(HzVal))+57.74462)


  elif freq == '600' and nuclei == 'N' and salt == "low":     
    # Construct: ??
    # Temperature: ??
    # Salt Conc: ??
    # Phosphate Conc: ??
    # Resonance: ??
    # Spectrometer: Bruker 600 (Duke)
    # Number of points: 44 ??
    # SL Cal Range: 25-3000 Hz ??
    # Error in duplicates: Yes
    # Date collected: ???
    # Originator: Yi X
    pldb23 = (-8.747*math.log(float(HzVal))+58.3631)
    
  elif freq == '700' and nuclei == 'N' and salt == "high":          
    # Construct: DickersonGT G/T.lbl
    # Temperature: 298K
    # pH: 5.93
    # Salt Conc: 125 mM NaCl
    # Phosphate Conc: 15 mM PO4
    # Resonance: Tdown-N3
    # Spectrometer: Bruker 700
    # Number of points: 25x1
    # SL Cal Range: 50-2000 Hz 
    # Error in duplicates: NA
    # Date collected: 03-26-2015
    # p1: 13.63
    # Fit type: Origin unweighted fit
    # Originator: Isaac K
    pldb23 = (-8.74862*math.log(float(HzVal))+57.93041)

    #Return the calculated power level
  return ('{0:.2f}'.format(round(pldb23, 2)))

## Old Calibrations Below ##
  # elif freq == '700' and nuclei == 'N' and salt == "low":         
  #   # Construct: wtTAR & sl-hpGT-GTG
  #   # Temperature: 298K
  #   # pH: 6.4 (TAR) and 8.4 (hpGT)
  #   # Salt Conc: 25 mM NaCl
  #   # Phosphate Conc: 15 mM
  #   # Resonance: TAR U37-N3 | hpGT G15-N1 & T5-N3
  #   #            Origin Weighted Fit
  #   # Spectrometer: Bruker 700
  #   # Number of points: 25 x3
  #   # SL Cal Range: 100-2000 Hz
  #   # Error in duplicates: Yes
  #   # Date collected: 01-11-2015
  #   # Originator: Isaac K
  #   pldb23 = (-8.73246*math.log(float(HzVal))+57.93763)

  # elif freq == '600' and nuclei == 'C' and salt == "low":
  #   # Construct: wtTAR
  #   # Temperature: 298K
  #   # pH: 6.4
  #   # Salt Conc: 25 mM NaCl
  #   # Phosphate Conc: 15 mM
  #   # Resonance: A35 C8?
  #   # Spectrometer: Bruker 600
  #   # Number of points: 50
  #   # SL Cal Range: 50-3500 Hz
  #   # Error in duplicates: None
  #   # Date collected: ???
  #   # Originator: Isaac K
  #   pldb23 = (-8.753*math.log(float(HzVal))+69.406)

  # if freq == '700' and nuclei == 'C' and salt == "low":
  #   # Construct: wtTAR 2
  #   # Temperature: 298K
  #   # pH: 6.4
  #   # Salt Conc: 25 mM NaCl
  #   # Phosphate Conc: 15 mM PO4
  #   # Resonance: A35 C8
  #   # Spectrometer: Bruker 700
  #   # Seq Type: 15N seq
  #   # Number of points: 25 x 2
  #   # SL Cal Range: 50-3500 Hz
  #   # Error in duplicates: Yes
  #   # Date collected: 01-18-2015
  #   # Originator: Isaac K / Dawn K
  #   pldb23 = (-8.74317*math.log(float(HzVal))+67.60758)

    # # Construct: wtTAR
    # # Temperature: 298K
    # # pH: 6.4
    # # Salt Conc: 25 mM NaCl
    # # Phosphate Conc: 15 mM
    # # Resonance: U41 N3 / U37 N3
    # #            Origin Weighted Fit
    # # Spectrometer: Bruker 700
    # # Number of points: 50 x2
    # # SL Cal Range: 100-2000 Hz
    # # Error in duplicates: Yes
    # # Date collected: ???
    # # Originator: Isaac K
    # pldb23 = (-8.73359*math.log(float(HzVal))+57.826)

  # elif freq == '700' and nuclei == 'C' and salt == "high":
  #   # Construct: cDGT
  #   # Temperature: 298K
  #   # pH: 6.8
  #   # Salt Conc: 125 mM NaCl
  #   # Phosphate Conc: 15 mM PO4
  #   # Resonance: T9 C6
  #   # Spectrometer: Bruker 700
  #   # Number of points: 50
  #   # SL Cal Range: 100-3500 Hz
  #   # Error in duplicates: None
  #   # Date collected: ???
  #   # Originator: Isaac K
  #   pldb23 = (-8.802*math.log(float(HzVal))+67.849)

  # elif freq == '700' and nuclei == 'N' and salt == "high":          
  #   # Construct: cDGT
  #   # Temperature: 298K
  #   # pH: 6.8
  #   # Salt Conc: 125 mM NaCl
  #   # Phosphate Conc: 15 mM PO4
  #   # Resonance: G4 N1
  #   # Spectrometer: Bruker 700
  #   # Number of points: 50
  #   # SL Cal Range: 100-2000 Hz
  #   # Error in duplicates: Yes
  #   # Date collected: ???
  #   # Originator: Isaac K
  #   pldb23 = (-8.727*math.log(float(HzVal))+57.79)

## 30 mM PO4 Buffer Calibration Below ##
  # elif freq == '700' and nuclei == 'N' and salt == "low":        
  #   # Construct: Dickerson G-T G/T-lbl
  #   # Temperature: 298K
  #   # pH: 7.02
  #   # Salt Conc: 0 mM NaCl
  #   # Phosphate Conc: 30 mM PO4
  #   # Resonance: G(downfield)-N1, T(upfield)-N3, T(downfield)-N3
  #   #            Origin Weighted Fit
  #   # Spectrometer: Bruker 700
  #   # Number of points: 25 x3
  #   # SL Cal Range: 100-2000 Hz
  #   # Error in duplicates: Yes
  #   # Date collected: Dec-14-2014
  #   # Originator: Isaac K
  #   pldb23 = (-8.73342*math.log(float(HzVal))+57.826)