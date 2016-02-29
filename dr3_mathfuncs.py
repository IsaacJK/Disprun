#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

# numpy imports
from numpy import absolute, array
from numpy import cos
from numpy import exp
from numpy import mean, std
from numpy import ones
from numpy.random import choice
from numpy.random import normal
import sys
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#--------------------------------------#
# Fit to single monoexponential decay
# # x = independent (time)
# # a = exponential prefactor
# # b = lambda (R1rho)
# # c = constant (ignore?)
# def func(x,a,b,c):
#   return a*exp(-b*x) + c
#--------------------------------------#
def expdec(x,a,b):
  return a*exp(-b*x)

#--------------------------------------#
# Fit to single damped sine decay
#  - R : relax in 1/sec
#  - omg : omega, angular freq in rad/sec
#  - phi : phase offset
#--------------------------------------#
def dampedsine(x,a,R,omg,phi):
  return(a*exp(-R*x)*cos(omg*x+phi))

#--------------------------------------#
# Calculate Chi-square for monoexponential decay
#--------------------------------------#
def redchi2(func, Params, xd, yd, yerr):
  chi2 = 0.

  chi2 += sum([((func(x,*Params)-y)/ye)**2. for x,y,ye in zip(xd,yd,yerr)])

  redchi2 = chi2 / (len(xd) - len(Params))
  
  return redchi2

#--------------------------------------#
# Calculate R1rho and standard error
#  from weighted L-M fit
#--------------------------------------#
def StandardError(func, p0, ints, delays, noiserel):

  # Calculate standard error from covariance matrix of L-M fit
  popt, pcov = curve_fit(func, delays, ints, p0=p0, sigma=noiserel)
  err = []

  # Square-root of diagonal elements of covariance matrix
  for i in range(len(popt)):
    try:
      err.append(absolute(pcov[i][i])**0.5)
    except:
      err.append(0.0)

  return popt, err

#--------------------------------------#
# MC Mean and Error in fitted function (2-param monoexp decay)
#  If StdErr = True, then Error is given by covariance matrix
#   of weighted L-M fit
#--------------------------------------#
def MonteCarlo(func, p0, ints, delays, noiserel, iterations, StdErr=False):
  iterations = int(iterations)
  oneList = ones(len(ints))
  # rndints = array([[normal(x,min(noiserel)) for x in ints] for y in range(iterations)])
  rndints = array([[normal(x,y*1.96) for x,y in zip(ints, noiserel)] for y in range(iterations)])

  # Fit random intensities MCnum of times
  # R1rhos = [curve_fit(func, delays, x, p0=p0)[0][1] for x in rndints]
  Fits = [curve_fit(func, delays, x, p0=p0, maxfev=25000) for x in rndints]

  # Split list of fits by A (prexponential factor), R (R1rho)
  #  and the covariance matrices
  FitVals, cov = [], []

  for i in Fits:
    FitVals.append(i[0])
    cov.append(i[1])
  # Cast them as arrays
  FitVals, cov = array(FitVals), array(cov)

  # Get MC-error in fitted values
  MCFit_mu, MCFit_err = [], []
  for v in FitVals.T:
    fit, err = norm.fit(v)
    MCFit_mu.append(fit)
    MCFit_err.append(err)

  # Calculate standard error from average covariance matrix
  if StdErr == True:
    meanCov = mean(cov, axis=0)
    err = []
    for i in range(len(p0)):
      try:
        err.append(absolute(meanCov[i][i])**0.5)
      except:
        err.append(0.)
    R1rho_std = err[1]

  return array(MCFit_mu), array(MCFit_err), array(FitVals)

#--------------------------------------#
# Random with replace Mean and Error in fitted function (2-param monoexp decay)
#  If StdErr = True, then Error is given by covariance matrix
#   of weighted L-M fit
#--------------------------------------#
def Bootstrap(func, p0, ints, delays, noiserel, iterations, StdErr=False):
  iterations = int(iterations)

  # Array of indices of length of ints array
  idxArr = [x for x in range(len(ints))]
  # Generate M number of idxArrays that have been randomized with replacement
  rndIdx = array([choice(idxArr, size=len(idxArr)) for _ in range(iterations)])
  # Generate array of delays,ints,noise for each random idxArray in rndIdx
  rndData = array([[delays[x],ints[x],noiserel[x]] for j in rndIdx for x in j])
  # Reshape to size MxNxO, M=iterations, N=len(ints), O=3 (delays, ints, noise)
  rndData = rndData.reshape((iterations,len(ints),3))

  Fits = [curve_fit(func, x[:,0], x[:,1], p0=p0, maxfev=10000) for x in rndData]
  
  # Split list of fits by A (prexponential factor), R (R1rho)
  #  and the covariance matrices
  FitVals, cov = [], []
  for i in Fits:
    FitVals.append(i[0])
    cov.append(i[1])
  # Cast them as arrays
  FitVals, cov = array(FitVals), array(cov)

  # Get MC-error in fitted values
  BSFit_mu, BSFit_err = [], []
  for v in FitVals.T:
    fit, err = norm.fit(v)
    BSFit_mu.append(fit)
    BSFit_err.append(err)

  # Calculate standard error from average covariance matrix
  if StdErr == True:
    meanCov = mean(cov, axis=0)
    err = []
    for i in range(len(p0)):
      try:
        err.append(absolute(meanCov[i][i])**0.5)
      except:
        err.append(0.)
    R1rho_std = err[1]

  return  array(BSFit_mu), array(BSFit_err), FitVals
