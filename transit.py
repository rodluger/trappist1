#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit.py
----------

The TRAPPIST-1 transit model class.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import pysyzygy as ps
import everest
import os
TRAPPIST_DAT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
G = 6.672e-8

def GetTransitTimes(file = 'ttv_kruse.dat'):
  '''

  '''

  planet, _, time, dtime = np.loadtxt(os.path.join(TRAPPIST_DAT, file), unpack = True)
  transit_times = [None for i in range(7)]
  if file == 'ttv_kruse.dat':
    for i in range(7):
      inds = np.where(planet == i + 1)[0]
      transit_times[i] = time[inds] + (2455000 - 2454833)  
  elif file == 'ttv_agol.dat':
    for i in range(6):
      inds = np.where(planet == i + 1)[0]
      transit_times[i] = time[inds] + (2450000 - 2454833)
    # Append a few extra for padding
    pad = [transit_times[i][-1] + np.median(np.diff(transit_times[i])),
           transit_times[i][-1] + 2 * np.median(np.diff(transit_times[i])),
           transit_times[i][-1] + 3 * np.median(np.diff(transit_times[i]))]
    transit_times[i] = np.append(transit_times[i], pad)

  return PlanetProperty(transit_times)
    
class PlanetProperty(object):
  '''
  
  '''
  
  def __init__(self, proplist):
    '''
    
    '''
    
    self.proplist = proplist
    self.planetnames = ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'][:len(self.proplist)]
    self.planetnumbers = list(range(len(self.planetnames)))
    self.propdict = {}
    for planet, item in zip(self.planetnames, self.proplist):
      self.propdict.update({planet: item})
  
  def __getitem__(self, key):
    '''
    
    '''
    
    if key in self.planetnames:
      return self.propdict[key]
    elif key in self.planetnumbers:
      return self.proplist[key]
    else:
      raise ValueError('Invalid value for `key`.')

  def __setitem__(self, key, value):
    '''
    
    '''
    
    if key in self.planetnames:
      self.propdict[key] = value
      self.proplist[np.argmax(np.array(self.planetnames) == key)] = value
    elif key in self.planetnumbers:
      self.proplist[key] = value
      self.propdict[self.planetnames[key]] = value
    else:
      raise ValueError('Invalid value for `key`.')

class Trappist1(object):
  '''
  
  '''
  
  def __init__(self):
    '''
    
    '''
    
    # Values from Gillon+ (2017)
    self.EPIC = 200164267
    self.planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    self.rhos = 50.7 * 1.4
    self.period = PlanetProperty([1.51087081, 2.4218233, 4.049610, 6.099615, 9.206690, 12.35294, 18.765])
    self.RpRs = PlanetProperty(np.sqrt([0.7266, 0.687, 0.367, 0.519, 0.673, 0.782, 0.352]) / 10.)
    self.sig_RpRs = PlanetProperty(np.sqrt([0.0088, 0.010, 0.017, 0.026, 0.023, 0.027, 0.0326]) / 10.)
    self.t0 = PlanetProperty(list(np.array([7322.51736, 7282.80728, 7670.141615, 7660.37859, 7671.39767, 7665.34937, 7662.55463]) + (2450000 - 2454833)))
    self.dur = PlanetProperty(np.array([36.40, 42.37, 49.13, 57.21, 62.60, 68.40, 76.7]) / 1440.)
    self.b = PlanetProperty([0.126, 0.161, 0.17, 0.12, 0.382, 0.421, 0.45])
    self.MpMs = PlanetProperty(list(np.array([0.85, 1.38, 0.41, 0.62, 0.68, 1.34, 0.]) * 3.74e-5))
    self.aRs = PlanetProperty([20.50, 28.08, 39.55, 51.97, 68.4, 83.2, 117.])
    self.inc = PlanetProperty([89.65, 89.67, 89.75, 89.86, 89.680, 89.710, 89.80])
    self.ecc = PlanetProperty([0., 0., 0., 0., 0., 0., 0.])
    self.w = PlanetProperty([90., 90., 90., 90., 90., 90., 90.])
    
    # TTVs
    self.times = GetTransitTimes()
    
    # Quadratic limb darkening table from Claret et al. (2011)
    # for filter = Kp, method = L, model = PHOENIX.
    # Data downloaded from `<CDS> http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/529/A75#sRM2.2`_.
    self.u1 = 1.0181
    self.u2 = -0.0404
      
  def transit_inds(self, time, planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h'], pad = 2.5):
    '''
    
    '''
    
    inds = []
    for planet in planets:
      
      # Get transit times
      times = self.times[planet]
      dur = self.dur[planet]
      period = self.period[planet]
      t0 = self.t0[planet]
      if times is None:
        t0 -= period * divmod(t0 - time[0], period)[0]
        n, r = divmod(time[0] - t0, period)
        if r < dur / 2.: t0 = t0 + n * period
        else: t0 = t0 + (n + 1) * period
        times = np.arange(t0, time[-1] + period, period)
      # Get indices
      for t in times:
        inds.extend(np.where(np.abs(time - t) <= (dur / 2) * pad)[0])

    return np.array(sorted(list(set(inds))), dtype = int)

  def flux(self, time, planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h'], cadence = 'lc'):
    '''
    
    '''
    
    # Ensure it's a list
    if type(planets) is str:
      planets = [planets]
    
    # Compute flux for each planet
    flux = np.ones_like(time)
    for planet in planets:
      if cadence == 'lc':
        model = ps.Transit(per = self.period[planet], b = self.b[planet], RpRs = self.RpRs[planet], t0 = self.t0[planet], 
                           rhos = self.rhos, ecc = self.ecc[planet], w = self.w[planet] * np.pi / 180., u1 = self.u1, 
                           u2 = self.u2, times = self.times[planet])
      else:
        model = ps.Transit(per = self.period[planet], b = self.b[planet], RpRs = self.RpRs[planet], t0 = self.t0[planet], 
                           rhos = self.rhos, ecc = self.ecc[planet], w = self.w[planet] * np.pi / 180., u1 = self.u1, 
                           u2 = self.u2, times = self.times[planet], exptime = ps.KEPSHRTCAD)
      flux *= model(time)
    
    return flux
  
  @property
  def transit_model(self):
    '''
    
    '''
    
    res = [i for i in range(7)]
    for i, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):
      res[i] = everest.TransitModel(planet, per = self.period[planet], b = self.b[planet], RpRs = self.RpRs[planet], t0 = self.t0[planet], 
                                    rhos = self.rhos, ecc = self.ecc[planet], w = self.w[planet] * np.pi / 180., u1 = self.u1, 
                                    u2 = self.u2, times = self.times[planet], sig_RpRs = self.sig_RpRs[planet])
    return res