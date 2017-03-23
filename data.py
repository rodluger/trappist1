#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
data.py
-------

The main interface to the K2 raw cadence data.

'''

from __future__ import division, print_function, absolute_import #, unicode_literals
from transit import Trappist1
import subprocess, os, shutil
import re
import numpy as np
import sys
import everest
from everest.utils import AP_COLLAPSED_PIXEL, AP_SATURATED_PIXEL
from everest.config import EVEREST_DAT
TRAPPIST_DAT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
TRAPPIST_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
try:
  import pyfits
except:
  import astropy.io.fits as pyfits
import matplotlib
from matplotlib.widgets import Slider
from matplotlib.ticker import FuncFormatter
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
import matplotlib.pyplot as pl
from scipy.ndimage import zoom
try:
  import pyfits
except ImportError:
  try:
    import astropy.io.fits as pyfits
  except ImportError:
    raise Exception('Please install the `pyfits` package.')
import logging
log = logging.getLogger(__name__)

# Some important EPIC IDs
TRAPPIST_EPIC = 200164267
SHORTCAD_NEIGHBORS = [245919787, 246011640, 246329409, 246331757, 246375295]
LONGCAD_NEIGHBORS = [246165150, 246211745, 246171759, 246127507, 246228828, 206392586,
                     246121678, 246229336, 246196866, 246217553, 246239441, 246144695]

# Some important constants
LONGCAD_BREAKPOINTS = [1173, 2320]

# The data directory
TRAPPIST_EVEREST_DAT = os.path.join(EVEREST_DAT, 'k2', 'c12', ('%09d' % TRAPPIST_EPIC)[:4] + '00000', ('%09d' % TRAPPIST_EPIC)[4:])

class Everest(everest.Everest):
  '''

  '''
  
  def __init__(self, fitsfile, ID = 0, quiet = False, clobber = False, cadence = 'lc', **kwargs):
    '''
    
    '''
    
    # Read kwargs
    self.ID = ID
    self._season = 12
    self.mission = 'k2'
    self.clobber = clobber
    self.cadence = cadence
    
    # Initialize preliminary logging
    if not quiet:
      screen_level = logging.DEBUG
    else:
      screen_level = logging.CRITICAL
    everest.utils.InitLog(None, logging.DEBUG, screen_level, False)

    # Load
    self.fitsfile = fitsfile
    self.model_name = pyfits.getheader(self.fitsfile, 1)['MODEL']
    self._weights = None
    self.load_fits()

class ApertureSelector(object):
  '''

  '''
  
  def __init__(self, time, images, title = 'Aperture'):
    '''
    
    '''
    
    self.cadence = 0
    self.time = time
    self.fig, self.ax = pl.subplots(1, figsize = (10,7))
    self.fig.subplots_adjust(left = 0.1, bottom = 0.25, top = 0.925, right = 0.45)
    self.images = images
    self.nt, self.ny, self.nx = self.images.shape
    self.x = np.arange(0, self.nx)
    self.y = np.arange(0, self.ny)
    self.aperture = np.zeros((self.ny, self.nx), dtype = int)
    self.aperture[self.ny // 2 - 2:self.ny // 2 + 2][:, self.nx // 2 - 2:self.nx // 2 + 2] = 1
    self.contour = None
    self.last_j = None
    self.last_i = None
    self.title = title
    
    # Slider
    self.axslider = pl.axes([0.105, 0.2, 0.34, 0.03])
    self.slider = Slider(self.axslider, '', 0, self.nt - 1, valinit=0, valfmt='%d')
    self.slider.valtext.set_x(0.5)
    self.slider.valtext.set_ha('center')
    self.slider.on_changed(self.replot)
    
    # Background
    self.axbkg = pl.axes([0.105, 0.05, 0.34, 0.125])
    bkg = self.colbkg
    self.bkgplot1, = self.axbkg.plot(self.x, bkg, 'ro')
    self.bkgplot2, = self.axbkg.plot(self.x, bkg, 'r-', alpha = 0.3)
    pad = 0.2 * (bkg.max() - bkg.min())
    self.axbkg.set_ylim(bkg.min() - pad, bkg.max() + pad)
    self.axbkg.set_xlim(-0.7, self.nx - 0.3)
    for tick in self.axbkg.get_yticklabels():
      tick.set_fontsize(7)
    self.axbkg.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p : '%.2f' % x))
    self.axbkg.set_ylabel('Bkg (%)', fontsize = 9)
    
    # Light curve
    self.axlc = pl.axes([0.5, 0.5, 0.4, 0.425])
    self.lcplot, = self.axlc.plot(self.time, self.flux, 'k.', alpha = 0.3, ms = 3)
    self.axlc.set_xticklabels([])
    self.axlc.yaxis.tick_right()
    self.axlc.set_ylabel('Light curve', fontsize = 14)
    self.lcstdtxt = self.axlc.annotate('%.2f ppm' % self.lcstd, xy = (0.025, 0.975), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 12, color = 'r')
    
    # Light curve background
    self.axlcbkg = pl.axes([0.5, 0.05, 0.4, 0.425])
    self.lcbkgplot, = self.axlcbkg.plot(self.time, self.lcbkg, 'k.', alpha = 0.3, ms = 3)
    self.axlcbkg.yaxis.tick_right()
    self.axlcbkg.set_ylabel('Background', fontsize = 14)
    self.bkgstdtxt = self.axlcbkg.annotate('%.2f ppm' % self.bkgstd, xy = (0.025, 0.975), xycoords = 'axes fraction', ha = 'left', va = 'top', fontsize = 12, color = 'r')
    
    # Trackers
    self.tracker1 = self.axlc.axvline(self.time[self.cadence], color = 'r', alpha = 0.5, lw = 1)
    self.tracker2 = self.axlcbkg.axvline(self.time[self.cadence], color = 'r', alpha = 0.5, lw = 1)
    
    # Appearance
    self.fig.canvas.set_window_title('Select an aperture')
    self.ax.axis('off')       
    self.ax.set_xlim(-0.7, self.nx - 0.3)
    self.ax.set_ylim(-0.7, self.ny - 0.3)
    self.ax.set_title(title, fontsize = 18)
    
    # Plot the image
    try:
      plasma = pl.get_cmap('plasma')
    except ValueError:
      plasma = pl.get_cmap('Greys')
    plasma.set_bad(alpha = 0)
    self.implot = self.ax.imshow(self.images[self.cadence], aspect = 'auto', interpolation = 'nearest', cmap = plasma, picker = True)
    self.fig.canvas.mpl_connect('motion_notify_event', self.mouse_drag)
    self.fig.canvas.mpl_connect('pick_event', self.mouse_click)
    
    # Update the contour
    self.update()
    
    # Enter interactive mode
    pl.show()
  
  @property
  def colbkg(self):
    '''
    
    '''
    
    # Flux in background pixels
    bkg = np.zeros(self.nx)
    for col in range(self.nx):
      b = np.where(self.aperture[:,col] == 0)
      bkg[col] = np.nanmedian(self.images[self.cadence][b,col])
    return 100 * (bkg / np.mean(bkg) - 1.)
  
  @property
  def lcbkg(self):
    '''
    
    '''
    
    binds = np.where(self.aperture ^ 1)
    bkg = np.nanmedian(np.array([f[binds] for f in self.images], dtype='float64'), axis = 1)
    return bkg.reshape(-1, 1)
  
  @property
  def flux(self):
    '''
    
    '''
    
    ap = np.where(self.aperture & 1)
    fpix2D = np.array([f[ap] for f in self.images], dtype='float64')
    return np.sum(fpix2D - self.lcbkg, axis = 1)
  
  @property
  def lcstd(self):
    '''
    
    '''
    
    return everest.k2.CDPP(self.flux)
  
  @property
  def bkgstd(self):
    '''
    
    '''
    
    return everest.k2.CDPP(self.lcbkg)
  
  def update_bkg(self):
    '''
    
    '''
    
    bkg = self.colbkg
    self.bkgplot1.set_ydata(bkg)
    self.bkgplot2.set_ydata(bkg)
    pad = 0.2 * (bkg.max() - bkg.min())
    self.axbkg.set_ylim(bkg.min() - pad, bkg.max() + pad)
    self.axbkg.set_xlim(-0.7, self.nx - 0.3)
  
  def update_lc(self):
    '''
    
    '''
    
    flux = self.flux
    self.lcplot.set_ydata(flux)
    pad = 0.2 * (flux.max() - flux.min())
    self.axlc.set_ylim(flux.min() - pad, flux.max() + pad)
    self.axlc.set_xlim(self.time[0], self.time[-1])
    self.lcstdtxt.set_text('%.2f ppm' % self.lcstd)
  
  def update_lcbkg(self):
    '''
    
    '''
    
    lcbkg = self.lcbkg
    self.lcbkgplot.set_ydata(lcbkg)
    pad = 0.2 * (lcbkg.max() - lcbkg.min())
    self.axlcbkg.set_ylim(lcbkg.min() - pad, lcbkg.max() + pad)
    self.axlcbkg.set_xlim(self.time[0], self.time[-1])
    self.bkgstdtxt.set_text('%.2f ppm' % self.bkgstd)
  
  def PadWithZeros(self, vector, pad_width, iaxis, kwargs):
    '''
    
    '''
    
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector
  
  def mouse_drag(self, event):
    '''
    
    '''
    
    if event.inaxes == self.ax and event.button == 1:

      # Index of nearest point
      i = np.nanargmin(((event.xdata - self.x) / self.nx) ** 2)
      j = np.nanargmin(((event.ydata - self.y) / self.ny) ** 2)  
      
      if (i == self.last_i) and (j == self.last_j):
        return
      else:
        self.last_i = i
        self.last_j = j
        
      # Toggle pixel
      if self.aperture[j,i]:
        self.aperture[j,i] = 0
      else:
        self.aperture[j,i] = 1
      
      # Update the contour
      self.update()

  def mouse_click(self, event):
    '''
    
    '''

    if event.mouseevent.inaxes == self.ax:
      
      # Index of nearest point
      i = np.nanargmin(((event.mouseevent.xdata - self.x) / self.nx) ** 2)
      j = np.nanargmin(((event.mouseevent.ydata - self.y) / self.ny) ** 2)  
      self.last_i = i
      self.last_j = j
        
      # Toggle pixel
      if self.aperture[j,i]:
        self.aperture[j,i] = 0
      else:
        self.aperture[j,i] = 1
      
      # Update the contour
      self.update()
    
  def update(self):
    '''
    
    '''
    
    # Update plot
    contour = np.zeros((self.ny,self.nx))
    contour[np.where(self.aperture)] = 1
    contour = np.lib.pad(contour, 1, self.PadWithZeros)
    highres = zoom(contour, 100, order = 0, mode='nearest') 
    extent = np.array([-1, self.nx, -1, self.ny])
    if self.contour is not None:
      for coll in self.contour.collections: 
        self.ax.collections.remove(coll) 
    self.contour = self.ax.contour(highres, levels=[0.5], extent=extent, origin='lower', colors='r', linewidths=2) 
    self.update_bkg()
    self.update_lc()
    self.update_lcbkg()
    self.fig.canvas.draw()

  def replot(self, val):
    '''
    
    '''
    
    # Update plot
    self.cadence = int(val)
    self.implot.set_data(self.images[int(val)])
    self.implot.set_clim(vmin = np.nanmin(self.images[int(val)]), vmax = np.nanmax(self.images[int(val)]))
    self.tracker1.set_xdata([self.time[self.cadence], self.time[self.cadence]])
    self.tracker2.set_xdata([self.time[self.cadence], self.time[self.cadence]])
    self.update_bkg()
    self.update_lc()
    self.update_lcbkg()
    self.fig.canvas.draw()

def GetData(EPIC, clobber_data = False, norm = 1.0, short_cadence = False, **kwargs):
  '''
  
  '''
  
  # Get the file name
  if short_cadence:
    filename = os.path.join(EVEREST_DAT, 'k2', 'c12', 
                           ('%09d' % EPIC)[:4] + '00000', ('%09d' % EPIC)[4:],
                           'data.sc.npz')
  else:
    filename = os.path.join(EVEREST_DAT, 'k2', 'c12', 
                           ('%09d' % EPIC)[:4] + '00000', ('%09d' % EPIC)[4:],
                           'data.npz')
  
  # Get tpf name
  if EPIC == TRAPPIST_EPIC:
    if short_cadence:
      tpf = os.path.join(TRAPPIST_DAT, 'trappist1-sc-tpf.fits.gz')
    else:
      tpf = os.path.join(TRAPPIST_DAT, 'trappist1-lc-tpf.fits.gz')
  else:
    tpf = os.path.join(TRAPPIST_DAT, 'ktwo%09d-unofficial-tpf.fits.gz' % EPIC)
  
  # Create the dir
  if not os.path.exists(os.path.dirname(filename)):
    os.makedirs(os.path.dirname(filename))
  
  # Check for saved data
  if not os.path.exists(filename) or clobber_data:

    print("Fetching data for target...")
    # Load the tpf
    f = pyfits.open(tpf)
    time = f[1].data['time']
    cadn = f[1].data['cadenceno']
    fpix = np.array(f[1].data['raw_cnts'], dtype = 'float64') / norm
    fpix[np.where(fpix < 0)] = np.nan
          
    # Get the static pixel images for plotting
    pixel_images = [fpix[0], fpix[len(fpix) // 2], fpix[len(fpix) - 1]]
  
    # Get the aperture interactively
    if short_cadence:
      aperture = ApertureSelector(time[::30], fpix[::30], title = 'EPIC %d' % EPIC).aperture
    else:
      aperture = ApertureSelector(time, fpix, title = 'EPIC %d' % EPIC).aperture
    if np.sum(aperture) == 0:
      raise ValueError("Empty aperture!")
    ap = np.where(aperture & 1)
    
    # Remove the background column by column
    bkg = np.zeros((fpix.shape[0], 1, fpix.shape[2]))      
    for col in range(fpix.shape[2]):
      binds = np.where(aperture[:,col] ^ 1)[0]
      bkg[:,0,col] = np.nanmedian(np.array([f[binds,col] for f in fpix], dtype='float64'), axis = 1)
    fpix -= bkg      
    bkg = np.nanmedian(bkg, axis = (1,2)).reshape(-1,1)

    # Saturation
    ncol = 0
    fpixnew = []
    saturation_tolerance = -0.1
    # This is an eyeballed value based on a handful of saturated stars...
    # Can only do better once we have the calibrated data!
    satflx = 9.e4 
    f97 = np.zeros((fpix.shape[1], fpix.shape[2]))
    for i in range(fpix.shape[1]):
      for j in range(fpix.shape[2]):
        if aperture[i,j]:
          # Let's remove NaNs...
          tmp = np.delete(fpix[:,i,j], np.where(np.isnan(fpix[:,i,j])))
          # ... and really bad outliers...
          if len(tmp):
            f = everest.math.SavGol(tmp)
            med = np.nanmedian(f)
            MAD = 1.4826 * np.nanmedian(np.abs(f - med))
            bad = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
            np.delete(tmp, bad)
            # ... so we can compute the 97.5th percentile flux
            i97 = int(0.975 * len(tmp))
            tmp = tmp[np.argsort(tmp)[i97]]
            f97[i,j] = tmp
    for j in range(aperture.shape[1]):
      if np.any(f97[:,j] > satflx):
        marked = False
        collapsed = np.zeros(len(fpix[:,0,0]))
        for i in range(aperture.shape[0]):
          if aperture[i,j]:
            if not marked:
              aperture[i,j] = AP_COLLAPSED_PIXEL
              marked = True
            else:
              aperture[i,j] = AP_SATURATED_PIXEL
            collapsed += fpix[:,i,j]
        if np.any(collapsed):
          fpixnew.append(collapsed)
          ncol += 1
      else:
        for i in range(aperture.shape[0]):
          if aperture[i,j]:
            fpixnew.append(fpix[:,i,j])
    fpix = np.array(fpixnew).T
    print("Collapsed %d saturated column(s)." % ncol)
    
    # SAP flux
    flux = np.sum(fpix, axis = 1)
  
    # Get NaN data points
    nanmask = np.where(np.isnan(flux) | (flux == 0) | np.isnan(time))[0]
  
    # Flag >10 sigma outliers
    t = np.delete(time, nanmask)
    f = np.delete(flux, nanmask)
    f = everest.math.SavGol(f)
    med = np.nanmedian(f)
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    bad = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
    badmask = np.array([np.argmax(time == t[i]) for i in bad], dtype = int)
  
    # Finalize the mask
    badmask = np.array(sorted(list(set(badmask))))
  
    # Interpolate the nans
    fpix = everest.math.Interpolate(time, nanmask, fpix)
    
    # Get metadata
    meta = [pyfits.getheader(tpf, 0).cards,
            pyfits.getheader(tpf, 1).cards,
            pyfits.getheader(tpf, 2).cards]
    
    # Save
    np.savez_compressed(filename, EPIC = EPIC, campaign = 12, time = time, fpix = fpix,
                        nanmask = nanmask, badmask = badmask, aperture = aperture, meta = meta, 
                        pixel_images = pixel_images, bkg = bkg, cadn = cadn)
  
  else:
    
    # Load
    foo = np.load(filename)
    time = foo['time']
    cadn = foo['cadn']
    fpix = foo['fpix']
    nanmask = foo['nanmask']
    badmask = foo['badmask']
    aperture = foo['aperture']
    meta = foo['meta']
    pixel_images = foo['pixel_images']
    bkg = foo['bkg']
  
  # Return
  data = everest.utils.DataContainer()
  data.ID = EPIC
  data.campaign = 12
  data.season = 12
  data.cadn = cadn
  data.time = time
  data.fpix = fpix
  data.fpix_err = np.sqrt(np.abs(fpix) * norm) / norm
  data.nanmask = nanmask
  data.badmask = badmask
  data.aperture = aperture
  data.aperture_name = 'manual'
  data.apertures = {'manual': aperture}
  data.quality = np.zeros_like(time, dtype = int)
  data.Xpos = None
  data.Ypos = None
  data.meta = meta
  if EPIC == TRAPPIST_EPIC:
    data.mag = 17.
  else:
    data.mag = np.nan
  data.pixel_images = pixel_images
  data.nearby = []
  data.hires = None
  data.saturated = False
  data.bkg = bkg
  return data

def Publish(star, cadence = 'lc'):
  '''
  
  '''
  
  everest.fits.MakeFITS(star)
  fitsfile = everest.missions.k2.FITSFile(star.ID, 12, cadence = cadence)
  filename = os.path.join(star.dir, star.name + '.fits')
  shutil.copy(os.path.join(star.dir, fitsfile), filename)
  os.remove(os.path.join(star.dir, fitsfile))
  return filename

def LongCadenceLightcurve(oiter = 20, osigma = 3, pad = 2.5, debug = True, mask_planets = True, clobber = False, **kwargs):
  '''
  The nPLD model with b-h transits masked prior to cross-validation.
  
  '''
  
  # Check if we already have a fits file
  fitsfile = os.path.join(TRAPPIST_EVEREST_DAT, 'nPLDTrappist.fits')
  if clobber or not os.path.exists(fitsfile):
  
    # Get neighbors
    neighbors_data = []
    for star in LONGCAD_NEIGHBORS:
      neighbors_data.append(GetData(star))
  
    # Get self
    data = GetData(TRAPPIST_EPIC)
  
    # Mask the transits
    if mask_planets:
      T1 = Trappist1()
      transitmask = T1.transit_inds(data.time, pad = pad)
    else:
      transitmask = np.array([], dtype = int)
    
    kwargs.update(dict(neighbors = LONGCAD_NEIGHBORS, 
                       neighbors_data = neighbors_data, 
                       breakpoints = LONGCAD_BREAKPOINTS,
                       oiter = oiter,
                       osigma = osigma,
                       debug = debug,
                       season = 12,
                       data = data))
  
    class nPLDTrappist(everest.nPLD):
      def load_tpf(self):
        super(nPLDTrappist, self).load_tpf()
        self.transitmask = transitmask
  
    # Run!
    model = nPLDTrappist(TRAPPIST_EPIC, **kwargs)
  
    # Generate the FITS file
    fitsfile = Publish(model)
  
  # Generate an emulated `Everest` instance
  star = Everest(fitsfile, ID = TRAPPIST_EPIC)
  
  return star

def ShortCadenceLightcurves():
  '''
  
  '''
  
  # Our Trappist-1 transit model instance
  T1 = Trappist1()
  times = T1.times['h']

  # Loop over the transits and generate a local SC de-trended light curve
  for i, t0 in enumerate(times):
  
    # Check if we've done this
    if os.path.exists(os.path.join(TRAPPIST_EVEREST_DAT, 'nPLDTrappisth%d.sc.fits' % (i + 1))):
      continue
    print("De-trending short cadence transit #%d..." % (i + 1))
    
    # Get SC data
    data = GetData(TRAPPIST_EPIC, short_cadence = True)
  
    # Apply a mask to get the light curve in the vicinity of the transit
    # For transit 3, let's get more of a baseline so we can fit out that flare!
    if i == 3:
      mask = np.where(np.abs(data.time - t0) < 1.)[0]
    elif i == 1:
      mask = np.where((data.time > 2941.3) & (data.time < 2942.5))[0]
    else:
      mask = np.where(np.abs(data.time - t0) < 0.75)[0]
  
    # Get neighbors
    neighbors_data = []
    for star in SHORTCAD_NEIGHBORS:
      ndata = GetData(star, short_cadence = True)
      
      # Paranoia check
      assert len(ndata.time) == len(data.time), "Time array mismatch!"
      
      # This is a bit hacky -- let's delete the data outside
      # the window we're de-trending in. This means we need
      # to re-caculate the indices of the outliers.
      ndata.time = ndata.time[mask]
      ndata.fpix = ndata.fpix[mask]
      ndata.fpix_err = ndata.fpix_err[mask]
      ndata.cadn = ndata.cadn[mask]
      f = everest.math.SavGol(np.sum(ndata.fpix, axis = 1))
      med = np.nanmedian(f)
      MAD = 1.4826 * np.nanmedian(np.abs(f - med))
      ndata.badmask = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
      ndata.nanmask = np.where(np.isnan(ndata.time))[0]
      neighbors_data.append(ndata)
  
    # Now, do the same thing for the TRAPPIST-1 data
    data.time = data.time[mask]
    data.fpix = data.fpix[mask]
    data.fpix_err = data.fpix_err[mask]
    data.cadn = data.cadn[mask]
    f = everest.math.SavGol(np.sum(data.fpix, axis = 1))
    med = np.nanmedian(f)
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    data.badmask = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
    data.nanmask = np.where(np.isnan(data.time))[0]
    
    # Mask transits of all the planets with a generous
    # 3.5x padding around them so we don't overfit anything
    transitmask = T1.transit_inds(data.time, pad = 3.5)
     
    if i == 1:
      transitmask = T1.transit_inds(data.time, pad = 3.5, planets = ['b', 'c', 'd', 'e', 'f', 'g'])
      transitmask = np.append(transitmask, np.where((data.time >= 2942.11) & (data.time < 2942.3))[0])
    
    # Mask the large instrumental(?) feature near the third transit
    if i == 2:
      data.badmask = np.append(data.badmask, np.where((data.time >= 2961.00) & (data.time < 2961.01))[0])
    
    # Mask some features near the fourth transit that could
    # screw up the PLD de-trending
    if i == 3:
      # Mask the weird feature
      data.badmask = np.append(data.badmask, np.where((data.time >= 2942.0) & (data.time < 2942.2))[0])
      # Mask the two flares and the transit of `h` -- note that this is a really wide mask,
      # so our de-trending power will be lower (but *unbiased*) during the transit of `h`
      data.badmask = np.append(data.badmask, np.where((data.time >= 2979.66) & (data.time < 2980.2))[0])
  
    # There's some pretty weird stuff going on in the light curve
    # in the vicinity of the large flare, and we get slightly better
    # results if we *don't* use the neighboring PLD regressors
    # in this region.
    if i == 3:
      neighbors = 0
      neighbors_data = []
    else:
      neighbors = SHORTCAD_NEIGHBORS
  
    # Define our custom `nPLD` instance
    class nPLDTrappist(everest.nPLD):
      def load_tpf(self):
        super(nPLDTrappist, self).load_tpf()
        self.transitmask = transitmask
  
    # That's it for hacks and fudges. Let's de-trend! Note the aggresive
    # outlier clipping and the large number of outlier iterations.
    star = nPLDTrappist(TRAPPIST_EPIC, season = 12, cadence = 'sc',
                        breakpoints = [], debug = True, 
                        data = data, oiter = 20, osigma = 3,
                        neighbors = neighbors, neighbors_data = neighbors_data,
                        pld_order = 3)

    # Make the FITS file
    everest.fits.MakeFITS(star)
  
    # Rename the files
    npzfile = os.path.join(star.dir, 'nPLDTrappist.sc.npz')
    pdffile = os.path.join(star.dir, 'nPLDTrappist.sc.pdf')
    fitsfile = everest.missions.k2.FITSFile(TRAPPIST_EPIC, 12, cadence = 'sc')
    os.rename(os.path.join(star.dir, fitsfile), os.path.join(star.dir, 'nPLDTrappisth%d.sc.fits' % (i + 1)))
    os.rename(os.path.join(star.dir, npzfile), os.path.join(star.dir, 'nPLDTrappisth%d.sc.npz' % (i + 1)))
    os.rename(os.path.join(star.dir, pdffile), os.path.join(star.dir, 'nPLDTrappisth%d.sc.pdf' % (i + 1)))
  
  # Load each of the light curves
  star = [None, None, None, None]
  for i, t0 in enumerate(times):
    star[i] = Everest(os.path.join(TRAPPIST_EVEREST_DAT, 'nPLDTrappisth%d.sc.fits' % (i + 1)), i + 1)
    # Delete the data we didn't de-trend
    sc_inds = np.where(star[i].time > 0)[0]
    star[i].time = star[i].time[sc_inds]
    star[i].fraw = star[i].fraw[sc_inds]
    star[i].fraw_err = star[i].fraw_err[sc_inds]
    star[i].cadn = star[i].cadn[sc_inds]
    star[i].model = star[i].model[sc_inds]
  
  return star