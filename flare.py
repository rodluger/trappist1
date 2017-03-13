#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 16:05:56 2017

@author: eakruse

A simple flare model based on Davenport et al. (2014, ApJ 797 122).

"""


import numpy as np
import matplotlib.pyplot as plt

def rise(thalf):
    return 1. + 1.941 * thalf - 0.175 * thalf**2. - 2.246 * thalf**3. - 1.125 * thalf**4.

def decay(thalf):
    return 0.689 * np.exp(-1.6 * thalf) + 0.303 * np.exp(-0.2783 * thalf)

def flare_model(times, tscale, peak, t0, a = 1, b = 0, c = 0, d = 0, e = 0):
    """
    times : ndarray
        Input array of times to calcluate the flare model

    tscale : float
        Characteristic time scale of the flare

    peak : float
        Maximum flux contributed by the flare. (0 for no flare)

    t0 : float
        Where the peak brightness of the flare occurs
        
    a-e : float
        Polynomial baseline parameters (optional)
        
    """
    outflux = np.ones(times.size)

    rtimes = times - t0

    if tscale < 0.:
      return np.ones_like(times)

    outflux[rtimes < -tscale] = 0.
    outflux[(rtimes >= -tscale) & (rtimes <= 0.)] = rise(rtimes[(rtimes >= -tscale) & (rtimes <= 0.)] / tscale)
    outflux[(rtimes >= 0.)] = decay(rtimes[(rtimes >= 0.)]/tscale)

    return (a + b * (times - t0) + c * (times - t0) ** 2 + d * (times - t0) ** 3 + e * (times - t0) ** 4) + np.array(outflux * peak)