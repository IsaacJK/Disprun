#########################################################################
# Disprun : R1rho RD Peak/Exp Fitting Program v3.06
#  Unfinished version.
#  Isaac Kimsey 09-17-2015
#
# Dependencies: numpy, matplotlib, scipy
#########################################################################

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import dr3_mathfuncs as mf
from numpy import absolute, linspace, array, mean, std, median
from scipy.stats import norm
import sys

def PlotDecay(dly, ints, noiserel, FitVals, R1p_mu, R1p_sigma, w1, wrf, fn, redchi2, p0, func):
  plt.close()
  fig, ax = plt.subplots(1,2, sharex=True, sharey=False, figsize=(16, 8))

  # Fine delays for trendline
  finedly = linspace(min(dly), max(dly), 50)

  # Label axes
  ax[0].set_xlabel("T-relax (sec)",fontsize=16)
  ax[1].set_xlabel("T-relax (sec)",fontsize=16)
  ax[0].set_ylabel("Intensity",fontsize=16)
  ax[1].set_ylabel("Residual Intensity",fontsize=16)
 
  # Set axes limits
  ax[0].set_xlim([min(dly)-max(dly)*0.05, max(dly)*1.05])
  ax[1].set_xlim([min(dly)-max(dly)*0.05, max(dly)*1.05])
  if min(ints) >= 0:
    ax[0].set_ylim([0., max(ints) * 1.05])
  else:
    ax[0].set_ylim([min(ints) - max(ints) * 0.05, max(ints) * 1.05])
  # Set axes titles
  ax[0].set_title(r'$R_{1\rho}\,\mathrm{%s\pm%s}\,s^{-1}\,|\,\omega_1\,%s\,Hz\,|\,\Omega_{eff}\,%s\,Hz\,|\,\overline{\chi}^2=%s$' 
                  % (round(R1p_mu,2),round(R1p_sigma,2), w1, wrf, round(redchi2,2)),
                  size=18)
  ax[1].set_title("Folder (%s)" % fn, size=18)

  tY, pY = [], []
  # Gather all trendlines
  for i in FitVals:
    tY = array([func(x, *i) for x in finedly])
    pY.append(tY)
    # ax[0].plot(finedly, func(finedly, *i),'-',c='red')

  # Plot best-fit, range of trendlines/barrier of error
  pY = array(pY)
  ax[0].plot(finedly, func(finedly, *p0), '-', c='red', linewidth=2)
  ax[0].fill_between(finedly, pY.max(axis=0), pY.min(axis=0), facecolor='deepskyblue', alpha=0.5)

  # Plot with Y-error
  ax[0].errorbar(dly, ints, fmt = 'o', yerr=noiserel, markersize=12, c='blue',
                 barsabove=True, ecolor='black')
  
  # Plot with Y-error
  ax[1].errorbar(dly, ints -func(dly, *p0), fmt = 'o', c='blue',
                 yerr=noiserel, markersize=12, barsabove=True, ecolor='black')

  fig.set_tight_layout(True)
  
  return fig

# Plot normal distribution of fitted R1rho values
# index corresponds to R1rho vals or other fit vals
def PlotDist(FitVals, R1p_mu, R1p_sigma, w1, wrf, fn, index=1):
  plt.close()
  fig = plt.figure()

  # Sort Fit values values
  R1rhos = FitVals[FitVals[:,index].argsort()][:,index]

  R1p_median = median(R1rhos)

  # Text for legend
  legStr = (r'$\mu=%.2f$\n$median=%.2f$\n$\sigma=%.2f$' 
            % (R1p_mu, R1p_median, R1p_sigma))
  props = dict(boxstyle='round')
  
  # Normal fit
  normFit = norm.pdf(R1rhos, mean(R1rhos), std(R1rhos))

  # Plot normal trend line
  plt.plot(R1rhos, normFit, "-", color="black")
  # Plot histogram
  plt.hist(R1rhos, normed=True, color="skyblue")

  # Set x-axes label
  plt.xlabel(r'$R_{1\rho}\,s^{-1}$', fontsize=16)

  # Set plot titles
  plt.suptitle(r'$\omega_1\,%s\,Hz\,|\,\Omega_{eff}\,%s\,Hz\,|\,FN\,%s$' % (w1, wrf, fn),
              size=16)

  plt.title(r'$\mu=%.2f\,|\,median=%.2f\,|\,\sigma=%.2f\,|\,n=%s$'
            % (R1p_mu, R1p_median, R1p_sigma, len(R1rhos)), size=14)


  # plt.text(0.05, 0.95, legStr, fontsize=14, verticalalignment='top', bbox=props)
  return fig