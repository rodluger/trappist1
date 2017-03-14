#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
trappist1.py
------------

Methods for plotting the TRAPPIST-1 K2 light curve and reproducing
several of the figures in the paper.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from transit import Trappist1
from data import LongCadenceLightcurve, ShortCadenceLightcurves, TRAPPIST_OUT, TRAPPIST_EVEREST_DAT
from flare import flare_model
import numpy as np
import shutil, os
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
try:
  from tqdm import tqdm
  prange = lambda x: tqdm(range(x))
except ImportError:
  prange = lambda x: range(x)
import logging
log = logging.getLogger(__name__)
BJDOFF = 2454833 - 2457700

def PlotFolded():
  '''
  Plots the long cadence light curve of TRAPPIST-1 folded
  on the transits of each of the seven planets.
  
  '''
  
  # Get the light curve
  star = LongCadenceLightcurve()
  
  # Define our transit model for all seven planets
  # We will do a joint instrumental/transit fit
  T1 = Trappist1()
  star.transit_model = T1.transit_model

  # There's a single (very) high outlier during one of the 
  # transits of h; let's mask it
  star.badmask = np.append(star.badmask, [3373])

  # Compute the joint model to get depths
  star.compute_joint()
  
  # Set up the plot
  fig = pl.figure(figsize = (5,8))
  axes = [pl.subplot2grid((4, 2), (0, 0), colspan = 1, rowspan = 1),
          pl.subplot2grid((4, 2), (0, 1), colspan = 1, rowspan = 1),
          pl.subplot2grid((4, 2), (1, 0), colspan = 1, rowspan = 1),
          pl.subplot2grid((4, 2), (1, 1), colspan = 1, rowspan = 1),
          pl.subplot2grid((4, 2), (2, 0), colspan = 1, rowspan = 1),
          pl.subplot2grid((4, 2), (2, 1), colspan = 1, rowspan = 1),
          pl.subplot2grid((4, 2), (3, 0), colspan = 2, rowspan = 1)]
  fig.subplots_adjust(hspace = 0.4)
  fig.suptitle('TRAPPIST-1', fontsize = 16)
  axes[-1].set_xlabel('Time from transit center (days)', fontsize = 14)
  
  # Plot each planet
  planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
  for planet, ax in zip(planets, axes):
    star.plot_transit_model(fold = planet, ax = ax, show = False) 
    ax.set_ylabel("")
    ax.set_xlabel("") 
    ax.set_title(planet, fontweight = 'bold')
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
      tick.set_fontsize(8)
  
  # Show!
  pl.show()

def DeltaChisq():
  '''
  Plot the delta chi-squared as a function of time for planet `h`. This 
  is Figure 5 in the paper.
  
  '''
  
  # Get the light curve
  star = LongCadenceLightcurve()
  
  # Define our transit model for all seven planets
  T1 = Trappist1()
  star.transit_model = T1.transit_model
  
  # Subtract the model fit for the first six planets
  star.fraw /= T1.flux(star.time, planets = ['b', 'c', 'd', 'e', 'f', 'g'])
  
  # Compute the delta-chisq
  time, depth, vardepth, delchisq = star.search(per = T1.period['h'], rhos = T1.rhos, u1 = T1.u1, u2 = T1.u2, b = T1.b['h'])
  delchisqcond = delchisq - (depth - 0.0035) ** 2 / vardepth
  fig, ax = pl.subplots(2, figsize = (10, 6))
  ax[0].plot(time + BJDOFF, delchisq, lw = 1, color = 'k', alpha = 0.65)
  ax[1].plot(time + BJDOFF, delchisqcond, lw = 1, color = 'k', alpha = 0.65)
  ax[0].set_ylabel(r'$\Delta \chi^2$', fontsize = 18, labelpad = 12)
  ax[1].set_ylabel(r'$\Delta \chi^2_\mathrm{cond}$', fontsize = 18)
  ax[1].set_xlabel('Barycentric Julian Date − 2,457,700 [day]', fontsize = 18)
  ax[0].set_xlim(time[0] + BJDOFF, time[-1] + BJDOFF)
  ax[1].set_xlim(time[0] + BJDOFF, time[-1] + BJDOFF)
  for n in range(5,9):
    t0 = T1.t0['h'] + n * T1.period['h']
    if n == 8:
      ttv = 0.025
    else:
      ttv = 0 
    dc = delchisq[np.argmin(np.abs(time-t0-ttv))]
    dcc = delchisqcond[np.argmin(np.abs(time-t0-ttv))]
    ax[0].annotate('', xy = (t0 + BJDOFF, dc), xytext = (t0 + BJDOFF, dc + 2.5), xycoords = 'data', 
                    textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'))
    ax[1].annotate('', xy = (t0 + BJDOFF, dcc), xytext = (t0 + BJDOFF, dcc + 2.5), xycoords = 'data', 
                    textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'))
  ax[0].annotate('FLARE', xy = (2920 + BJDOFF, 6), ha = 'center', va = 'center', color = 'k', alpha = 0.5)
  ax[0].annotate('FLARE', xy = (2926.92 + BJDOFF, 6), ha = 'center', va = 'center', color = 'k', bbox=dict(fc='w',ec='w',alpha=0.75), alpha = 0.5)
  ax[0].annotate('FLARE', xy = (2947.16 + BJDOFF, 3.65), ha = 'center', va = 'center', color = 'k', bbox=dict(fc='w',ec='w',alpha=0.75), alpha = 0.5)
  ax[0].annotate('FLARE', xy = (2980.3 + BJDOFF, 6.3), ha = 'center', va = 'center', color = 'k', alpha = 0.5)
  ax[0].annotate('FLARE', xy = (2973.68 + BJDOFF, 3.65), ha = 'center', va = 'center', color = 'k', alpha = 0.5)
  ax[1].annotate('FLARE', xy = (2947.16 + BJDOFF, -5), ha = 'center', va = 'center', color = 'k', alpha = 0.5)
  ax[1].annotate('FLARE', xy = (2920.13 + BJDOFF, -5), ha = 'center', va = 'center', color = 'k', bbox=dict(fc='w',ec='w',alpha=0.75), alpha = 0.5)
  ax[1].annotate('FLARE', xy = (2926.9 + BJDOFF, -5), ha = 'center', va = 'center', color = 'k', alpha = 0.5)
  ax[1].annotate('FLARE', xy = (2980.3 + BJDOFF, -5), ha = 'center', va = 'center', color = 'k', bbox=dict(fc='w',ec='w',alpha=0.75), alpha = 0.5)
  ax[1].annotate('FLARE', xy = (2973.68 + BJDOFF, -5), ha = 'center', va = 'center', color = 'k', bbox=dict(fc='w',ec='w',alpha=0.75), alpha = 0.5)
  
  # Show
  pl.show()
   
def PowerSpectrum():
  '''
  Plot the delta chi-squared power spectrum as a function of period for planet `h`. This 
  is Figure 6 in the paper.
  
  '''
  
  # Get the light curve
  star = LongCadenceLightcurve()
  
  # Define our transit model for all seven planets
  T1 = Trappist1()
  star.transit_model = T1.transit_model
  
  # Subtract the model fit for the first six planets
  star.fraw /= T1.flux(star.time, planets = ['b', 'c', 'd', 'e', 'f', 'g'])
  
  # Compute the delta-chisq
  time, depth, vardepth, delchisq = star.search(per = T1.period['h'], rhos = T1.rhos, u1 = T1.u1, u2 = T1.u2, b = T1.b['h'])
  delchisqcond = delchisq - (depth - 0.0035) ** 2 / vardepth
  
  # Loop over a grid of periods
  periods = np.linspace(1, 50, 500000)
  P = np.zeros_like(periods)
  Pttvs = np.zeros_like(periods)
  lt = np.zeros_like(periods)
  t0 = T1.t0['h']
  dur = T1.dur['h']
  for i in prange(len(periods)):
    period = periods[i]
    t0i = t0 - period * divmod(t0 - time[0], period)[0]
    n, r = divmod(time[0] - t0i, period)
    if r < dur / 2.: t0i = t0i + n * period
    else: t0i = t0i + (n + 1) * period
    times = np.arange(t0i, time[-1], period)
    lt[i] = times[-1]
  
    # No TTVS
    P[i] = np.nansum([np.interp(t, time, delchisqcond) for t in times])
  
    # With TTVS
    for t in times:
      ind = np.argmin(np.abs(time - t)) 
      if ind > 1: 
        foo = np.nanmax(delchisqcond[ind-1:ind+1])
        if np.isfinite(foo):
          Pttvs[i] += foo

  # Ensure the grid is fine enough
  assert np.max(np.diff(lt)) < 0.02, "Period grid is not fine enough!"

  # Plot
  fig, ax = pl.subplots(2, sharex = True, figsize = (8, 5))
  ax[0].plot(periods, P, color = 'k', alpha = 0.75, lw = 1, zorder = -1)
  ax[1].plot(periods, Pttvs, color = 'k', alpha = 0.75, lw = 1, zorder = -1)
  ax[0].set_rasterization_zorder(0)
  ax[1].set_rasterization_zorder(0)
  ax[0].set_ylabel(r'$\Delta \chi^2_\mathrm{P}$ (periodic)', fontsize = 16)
  ax[1].set_ylabel(r'$\Delta \chi^2_\mathrm{P}$ (TTVs)', fontsize = 16)
  ax[1].set_xlabel(r'Period [days]', fontsize = 16)
  ax[0].set_ylim(-65,20)
  ax[1].set_ylim(-65,20)
  ax[0].annotate('h', xy = (18.765, -12), xytext = (18.765, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[0].annotate(r'$\frac{1}{2}$h', xy = (18.765 / 2, -18), xytext = (18.765 / 2, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[0].annotate(r'2h', xy = (18.765 * 2, -9), xytext = (18.765 * 2, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[0].annotate(r'$\frac{1}{3}$h', xy = (18.765 / 3, -26), xytext = (18.765 / 3, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[1].annotate('h', xy = (18.765, -10), xytext = (18.765, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[1].annotate(r'$\frac{1}{2}$h', xy = (18.765 / 2, -17), xytext = (18.765 / 2, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[1].annotate(r'$\frac{1}{3}$h', xy = (18.765 / 3, -22), xytext = (18.765 / 3, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  ax[1].annotate(r'2h', xy = (18.765 * 2, -7), xytext = (18.765 * 2, -50), xycoords = 'data', 
                 textcoords = 'data', arrowprops=dict(arrowstyle="->",color='r'),
                 ha = 'center', va = 'center', color = 'r')
  
  # Show
  pl.show()

def ShortCadence():
  '''
  Plot the short cadence data and diagnostics plots for
  TRAPPIST-1h. These are figures 2, 3, and 4 in the paper.
    
  '''
  
  # A Trappist-1 transit model instance
  T1 = Trappist1()

  # Load the long cadence light curve
  lc_star = LongCadenceLightcurve()
  
  # Load the short cadence light curves
  sc_stars = ShortCadenceLightcurves()
  
  # Set up our figure
  fig, ax = pl.subplots(3,2, figsize = (8, 8))
  fig.subplots_adjust(wspace = 0.15, hspace = 0.15)
  ax = ax.flatten()

  # The folded time arrays
  sc_time_fold = []
  sc_flux_fold = []
  lc_time_fold = []
  lc_flux_fold = []

  # Loop over the panels in the figure
  for i, n, label, t0 in zip(range(6), [1,2,3,4,3,4], ['1','2','3a','4a','3b','4b'], list(T1.times['h']) + list(T1.times['h'])[2:]):

    # Get this portion of the SC light curve
    sc_star = sc_stars[n - 1]
    sc_time = sc_star.time
    sc_flux = sc_star.flux / np.nanmedian(sc_star.flux)

    # Get this portion of the LC light curve
    lc_inds = np.where(np.abs(lc_star.time - t0) < 0.75)[0]
    lc_time = lc_star.time[lc_inds]
    lc_flux = lc_star.flux[lc_inds] / np.nanmedian(lc_star.flux[lc_inds])
  
    # HACK: Fit out the simultaneous transit of `b`
    if i == 4:
      tmp = np.array(T1.times['b'])
      T1.times['b'] = [2960.91]
      lc_transit_model = T1.flux(lc_time, planets = ['b'])
      sc_transit_model = T1.flux(sc_time, planets = ['b'], cadence = 'sc')
      
      # Create a little figure with the double transit model
      figbh, axbh = pl.subplots(3, figsize = (5,8))
      figbh.subplots_adjust(left = 0.2, top = 0.93)
      figbh.suptitle('Transit #3', fontsize = 16)
      for j in range(3):
        axbh[j].plot(sc_time + BJDOFF, sc_flux, 'k.', alpha = 0.75, ms = 3)
        axbh[j].set_xlim(t0 - 0.08 + BJDOFF, t0 + 0.08 + BJDOFF)
        axbh[j].set_ylim(0.965, 1.025)
        axbh[j].set_ylabel('Normalized Flux', fontsize = 12)
        axbh[j].get_xaxis().set_major_locator(MaxNLocator(4))
      axbh[0].set_xticklabels('')
      axbh[1].set_xticklabels('')
      axbh[2].set_xlabel('Barycentric Julian Date − 2,457,700 [day]', fontsize = 12)
      thires = np.linspace(t0 - 0.08, t0 + 0.08, 1000)
      axbh[0].plot(thires + BJDOFF, T1.flux(thires, planets = ['b']), ls = '-', color = 'orange')
      axbh[1].plot(thires + BJDOFF, T1.flux(thires, planets = ['h']), ls = '-', color = 'blue')
      axbh[2].plot(thires + BJDOFF, T1.flux(thires, planets = ['b', 'h']), 'r-')
      axbh[0].annotate('b', xy = (93.91, 1.01), xycoords = 'data', ha = 'center', va = 'center', fontsize = 18, color = 'orange')
      axbh[1].annotate('h', xy = (t0 + BJDOFF, 1.01), xycoords = 'data', ha = 'center', va = 'center', fontsize = 18, color = 'blue')
      axbh[2].annotate('b + h', xy = (93.915, 1.01), xycoords = 'data', ha = 'center', va = 'center', fontsize = 18, color = 'r')
      
      # Divide out the model for `b`
      lc_flux = lc_flux / lc_transit_model
      sc_flux = sc_flux / sc_transit_model
      T1.times['b'] = np.array(tmp)

    # HACK: Fit out the flare during the last transit
    if i == 5:
      x = np.array(sc_time)
      y = np.array(sc_flux)
      ii = np.where((x - t0 < 0.15) & (x - t0 > -0.275))[0]
      x = x[ii]
      y = y[ii]
      y[np.where(y < 0.98)] = 1.
      popt, pcov = curve_fit(flare_model, x, y, p0 = (0.005, 0.07, 2979.66, 1))
    
      # Create a little figure with the flare fit
      figflare, axflare = pl.subplots(1)
      figflare.suptitle('Transit #4', fontsize = 16)
      axflare.plot(sc_time + BJDOFF, sc_flux, 'k.', alpha = 0.75, ms = 3)
      thires = np.linspace(t0 - 0.275, t0 + 0.15, 1000)
      axflare.plot(thires + BJDOFF, flare_model(thires,*popt), 'r-', alpha = 0.3, ms = 3)
      axflare.set_xlim(-0.15 + t0 + BJDOFF, 0.15 + t0 + BJDOFF)
      axflare.set_ylim(0.97, 1.1)
      axflare.set_xlabel('Barycentric Julian Date − 2,457,700 [day]', fontsize = 14)
      axflare.annotate('h', xy = (t0 + BJDOFF, 1.01), xycoords = 'data', ha = 'center', va = 'center', fontsize = 18, color = 'r')
      axflare.set_ylabel('Normalized Flux', fontsize = 14)
      
      # Correct the arrays
      sc_flux /= flare_model(sc_time, *popt)
      lc_time_hires = np.linspace(lc_time[0],lc_time[-1],len(lc_time)*30)
      lc_flare_hires = flare_model(lc_time_hires, *popt)
      lc_flare = np.mean(lc_flare_hires.reshape(-1, 30), axis=1)
      lc_flux /= lc_flare
  
    # Folded arrays
    if i in [0,1,4,5]:
      sc_time_fold.extend(sc_time - t0)
      sc_flux_fold.extend(sc_flux)
      lc_time_fold.extend(lc_time - t0)
      lc_flux_fold.extend(lc_flux)
          
    # Bin SC to LC
    pad = 30 - (len(sc_time) % 30)
    if pad < 10:
      pad += 30
    offset = 10
    sc_time_bin = np.mean(np.concatenate([np.zeros(offset), sc_time, np.zeros(pad - offset)]).reshape(-1, 30), axis=1)
    sc_flux_bin = np.mean(np.concatenate([np.zeros(offset), sc_flux, np.zeros(pad - offset)]).reshape(-1, 30), axis=1)
    inds = np.where(np.abs(sc_time_bin - t0) < 0.275)
    sc_flux_bin /= np.nanmedian(sc_flux_bin[inds])
    sc_time_bin = sc_time_bin[1:-1]
    sc_flux_bin = sc_flux_bin[1:-1]
  
    # There's some nasty outliers in the LC data to the right of
    # transit #4; let's mask them to reduce clutter in the plot
    if i == 3 or i == 5:
      lc_flux[np.where(lc_time - t0 > 0.161)] = np.nan
      sc_flux_bin[np.where(sc_time_bin - t0 > 0.161)] = np.nan
  
    # Plot!
    ax[i].plot(sc_time - t0, sc_flux, 'k.', alpha = 0.3, ms = 3)
    ax[i].plot(sc_time_bin - t0, sc_flux_bin, color = 'orange', ls = '-', alpha = 1)
    ax[i].plot(lc_time - t0, lc_flux, 'r-', alpha = 0.75)
    
    # Lims and labels
    ax[i].set_ylim(0.97, 1.03)
    ax[i].set_xlim(-0.275, 0.275)
    if i == 4 or i == 5:
      ax[i].set_xlabel('Time from transit center (days)', fontsize = 13)
    if i == 0 or i == 2 or i == 4:
      ax[i].set_ylabel('Normalized flux', fontsize = 13)
    ax[i].annotate(label, xy = (0.05, 0.95), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 10, fontweight = 'bold')

  # Mark the other visible planets
  ax[0].annotate('b', xy = (-0.25, 1.01), color = 'r', ha = 'center')
  ax[2].annotate('c', xy = (-0.11, 1.01), color = 'r', ha = 'center')
  ax[2].annotate('b,h', xy = (0, 1.01), color = 'r', ha = 'center')
  ax[3].annotate('e', xy = (0.155, 1.01), color = 'r', ha = 'center')
  ax[4].annotate('c', xy = (-0.1, 1.01), color = 'r', ha = 'center')
  ax[5].annotate('e', xy = (0.155, 1.01), color = 'r', ha = 'center')

  # Mark our guesses
  ax[0].annotate('h', xy = (0, 1.01), color = 'r', ha = 'center')
  ax[1].annotate('h', xy = (0., 1.01), color = 'r', ha = 'center')
  ax[3].annotate('h', xy = (0., 1.01), color = 'r', ha = 'center')
  ax[4].annotate('h', xy = (0., 1.01), color = 'r', ha = 'center')
  ax[5].annotate('h', xy = (0., 1.01), color = 'r', ha = 'center')

  # Plot folded
  figfold = pl.figure(figsize = (7,5))
  axf = pl.subplot2grid((3, 1), (0, 0), colspan = 1, rowspan = 2)
  axr = pl.subplot2grid((3, 1), (2, 0), colspan = 1, rowspan = 1)
  
  # Short cadence
  sc_time_fold = np.array(sc_time_fold)
  sc_flux_fold = np.array(sc_flux_fold)
  inds = np.argsort(sc_time_fold)
  sc_time_fold = sc_time_fold[inds]
  sc_flux_fold = sc_flux_fold[inds]
  pad = 30 - (len(sc_time_fold) % 30)
  offset = 10
  sc_flux_fold_clipped = np.array(sc_flux_fold)
  sc_flux_fold_clipped[np.where(sc_flux_fold > 1.015)] = 1.
  sc_flux_fold_clipped[np.where(sc_flux_fold < 0.98)] = 1.
  sc_time_bin = np.mean(np.concatenate([np.zeros(offset), sc_time_fold, np.zeros(pad - offset)]).reshape(-1, 30), axis=1)
  sc_flux_bin = np.mean(np.concatenate([np.zeros(offset), sc_flux_fold_clipped, np.zeros(pad - offset)]).reshape(-1, 30), axis=1)
  sc_flux_bin /= np.nanmedian(sc_flux_bin)
  sc_time_bin = sc_time_bin[1:-1]
  sc_flux_bin = sc_flux_bin[1:-1]

  # Long cadence
  lc_time_fold = np.array(lc_time_fold)
  lc_flux_fold = np.array(lc_flux_fold)
  inds = np.argsort(lc_time_fold)
  lc_time_fold = lc_time_fold[inds]
  lc_flux_fold = lc_flux_fold[inds]
  lc_flux_fold[np.where(lc_flux_fold > 1.015)] = np.nan
  lc_flux_fold[np.where(lc_flux_fold < 0.98)] = np.nan
  
  # Plot it
  x = np.linspace(-0.35, 0.35, 1000)
  y = T1.flux(x + T1.times['h'][0], planets='h', cadence='sc')
  axf.plot(sc_time_fold,sc_flux_fold,'k.',alpha = 0.35, ms = 3, label = 'Short cadence')
  axf.plot(sc_time_bin,sc_flux_bin,color = 'orange', ls = '-', alpha = 1, label = 'Binned', lw = 2)
  axf.plot(x,y,'r-', alpha = 0.75, label = 'Transit model', lw = 2)
  axf.set_xlim(-0.35, 0.35)
  axf.set_ylim(0.96, 1.04)
  
  # Plot the transit model and the residuals
  sc_model = T1.flux(sc_time_fold + T1.times['h'][0], planets='h', cadence='sc')
  axr.plot(sc_time_fold, sc_flux_fold - sc_model, 'k.',alpha = 0.35, ms = 3)
  sc_flux_fold_clipped = np.array(sc_flux_fold - sc_model)
  sc_flux_fold_clipped[np.where(sc_flux_fold - sc_model > 0.015)] = 0.
  sc_flux_fold_clipped[np.where(sc_flux_fold - sc_model < -0.02)] = 0.
  sc_time_binr = np.mean(np.concatenate([np.zeros(offset), sc_time_fold, np.zeros(pad - offset)]).reshape(-1, 30), axis=1)
  sc_flux_binr = np.mean(np.concatenate([np.zeros(offset), sc_flux_fold_clipped, np.zeros(pad - offset)]).reshape(-1, 30), axis=1)
  sc_flux_binr -= np.nanmedian(sc_flux_binr)
  sc_time_binr = sc_time_binr[1:-1]
  sc_flux_binr = sc_flux_binr[1:-1]
  axr.plot(sc_time_binr,sc_flux_binr,color = 'orange', ls = '-', alpha = 1, lw = 2)
  axr.set_ylim(-0.02, 0.02)
  axr.set_xlim(-0.35, 0.35)
  
  # Final tweaks!
  axf.legend(loc = 'upper left', numpoints=3)
  axf.set_ylabel('Folded Normalized Flux', fontsize = 14, labelpad = 10)
  axf.set_xticklabels([])
  axr.set_ylabel('Residuals', fontsize = 14)
  axr.set_xlabel('Time from transit center (days)', fontsize = 14)
  
  # Show all plots
  pl.show()

if __name__ == '__main__':
  
  # Create the TRAPPIST1 folder
  if not os.path.exists(TRAPPIST_EVEREST_DAT):
    os.makedirs(TRAPPIST_EVEREST_DAT)
  fitsfiles = ['nPLDTrappist.fits', 'nPLDTrappisth1.sc.fits', 'nPLDTrappisth2.sc.fits', 
               'nPLDTrappisth3.sc.fits', 'nPLDTrappisth4.sc.fits']
  
  # NOTE: Comment these lines if you want to do the de-trending yourself!
  for file in fitsfiles:
    if not os.path.exists(file):
      shutil.copy(os.path.join(TRAPPIST_OUT, file), TRAPPIST_EVEREST_DAT)
  
  # Folded light curve for all planets
  PlotFolded()
  
  # Figures 2, 3, 4
  ShortCadence()
  
  # Figure 5
  DeltaChisq()
  
  # Figure 6
  PowerSpectrum()