<p align="center">
  <img width = "300" src="http://staff.washington.edu/rodluger/everest/_images/everest.png"/>
</p>
<p align="center">
  <a href="#"><img src="https://img.shields.io/badge/arXiv-XXXX.YYYYY-blue.svg?style=flat"/></a>
  <a href="https://raw.githubusercontent.com/rodluger/trappist1/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-brightgreen.svg"/></a>
  <a href="https://doi.org/10.5281/zenodo.376863"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.376863.svg"/></a>
</p>


## Installation
In order to run the scripts and interact with the light curves in this repository, you will need to install the latest <b>EVEREST</b> release (>=2.0.7):
<br/>
<pre><code>pip install everest-pipeline</code></pre>
If you alread have <b>EVEREST</b> installed, please upgrade it:
<br/>
<pre><code>pip install everest-pipeline --upgrade</code></pre>
<br/>
For more information on installing and using <b>EVEREST</b>, check out the [project github page](https://github.com/rodluger/everest). 

The methods in [trappist1.py](trappist1.py) allow users to plot and interact with the <b>EVEREST</b> light curve for TRAPPIST-1, as well as to reproduce several of the figures in the <a href="#">paper</a>.

## trappist1.PlotFolded()
Folded long cadence plots for each of the seven planets transiting TRAPPIST-1:
<p align="center">
  <img src="output/folded.png" width="50%"/>
</p>

## trappist1.ShortCadence()
The long cadence data folded on the period of planet <b>h</b>:
<p align="center">
  <img src="output/sc_folded.png" width="50%"/>
</p>
Each of the individual four transits of planet <b>h</b>, with a simultaneous transit of <b>b</b> and a near-simultaneous flare removed in the bottom panels:
<p align="center">
  <img src="output/sc_transits.png" width="50%"/>
</p>
A closer look at what's going on during the third transit:
<p align="center">
  <img src="output/transit_3.png" width="50%"/>
</p>
A closer look at what's going on during the fourth transit, with a flare fit based on [Davenport et al. (2014)](http://adsabs.harvard.edu/abs/2014ApJ...797..122D):
<p align="center">
  <img src="output/transit_4.png" width="50%"/>
</p>

## trappist1.DeltaChisq()
The delta-chi squared long cadence plot (top) and the delta-chi squared conditioned on the true depth of planet <b>h</b> (bottom):
<p align="center">
  <img src="output/deltachisq.png" width="50%"/>
</p>

## trappist1.PowerSpectrum()
The delta-chi squared long cadence power spectrum, where the period of <b>h</b> and its aliases are clearly visible:
<p align="center">
  <img src="output/powerspec.png" width="50%"/>
</p>
