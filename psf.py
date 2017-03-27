#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
psf.py
------

We use a 6x6 square aperture centered on TRAPPIST-1 to detrend the K2
pixel light curve with EVEREST. We use this script to calculate the
approximate amount of flux loss from this aperture by fitting a simple
integrated Gaussian (i.e., an error function) to the stellar PSF at
different points throughout the timeseries. We find that at times
when the star is closest to the edge of the aperture, the flux loss
is on the order of 4e-4, which for our purposes is negligible.

'''

from __future__ import division, print_function, absolute_import #, unicode_literals
from data import Everest, TRAPPIST_EPIC
import numpy as np
import matplotlib.pyplot as pl
from scipy.special import erf
from scipy.optimize import fmin
try:
  from tqdm import tqdm
  prange = lambda x: tqdm(range(x))
except ImportError:
  prange = lambda x: range(x)
  
def Erf2D(x, y, xc, yc, a, b, sigma):
  '''
  A simple 2D error function PSF model.
  
  '''
  c = np.sqrt(2) * sigma
  ex = (erf((x + 1 - xc) / c) - erf((x - xc) / c))
  ey = (erf((y + 1 - yc) / c) - erf((y - yc) / c))
  return a * ex * ey + b

def ChiSq(params, fpix):
  '''
  Erf fit residuals for optimization.
  
  '''
  
  xc, yc, a, b, sigma = params
  y, x = np.meshgrid(np.arange(6),np.arange(6))
  return np.sum((Erf2D(x,y,xc,yc,a,b,sigma) - fpix) ** 2)

def GetSigma():
  '''
  Get the standard deviation of the Gaussian PSF. I find `sigma` = 0.45 pixels.
  I also plotted the `x` and `y` obtained below to find that average the maximum extent of
  the stellar motion is x,y ~ y,x ~ (2.5, 3.5).
  
  '''
  
  star = Everest('/Users/rodrigo/src/trappist1/output/nPLDTrappist.fits', TRAPPIST_EPIC)
  fpix = star.fpix.reshape(-1, 6, 6).swapaxes(1,2)
  guess = [3., 3., 1e4, 1e2, 0.5]
  n = 0
  niter = len(fpix) // 10
  x = np.zeros(niter)
  y = np.zeros(niter)
  a = np.zeros(niter)
  b = np.zeros(niter)
  sigma = np.zeros(niter)
  for n in prange(niter):
    x[n], y[n], a[n], b[n], sigma[n] = fmin(ChiSq, guess, args = (fpix[n * 10],), disp = 0)
  return np.nanmedian(sigma)

# Aperture center
x, y = np.meshgrid(np.arange(3 - 10, 3 + 10),np.arange(3 - 10, 3 + 10))
image = Erf2D(x, y, 3, 3, 1, 0, 0.45)
aperture = np.where((x > 0) & (x < 6) & (y > 0) & (y < 6))
influx = np.sum(image[aperture])
outflux = np.sum(image) - influx
print("Flux loss @ aperture center: %.2e" % (outflux / (influx + outflux)))

# Maximum deviation from aperture center
x, y = np.meshgrid(np.arange(3 - 10, 3 + 10),np.arange(3 - 10, 3 + 10))
image = Erf2D(x, y, 2.5, 3.5, 1, 0, 0.45)
aperture = np.where((x > 0) & (x < 6) & (y > 0) & (y < 6))
influx = np.sum(image[aperture])
outflux = np.sum(image) - influx
print("Flux loss @ maximum deviation: %.2e" % (outflux / (influx + outflux)))