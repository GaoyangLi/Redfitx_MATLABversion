# Redfix_MATLABversion
coherency for unevenly sampled time series
redfit_x.pdf, tauest.pdf: some theoretical background;
also refer to http://adsbit.harvard.edu/full/1989ApJ...343..874S
for Scargle-Lomb spectrum;

brent.m: Brent method for locating a local minimum;
ls.m: the function to be minimised in tau estimation;
redfitx.m: (cross) spectral estimation of two unevenly sampled datasets
segmenting.m: function that cuts time-series into segments for windowed
              spectral analysis;
SLspectrum.m: Computing Scargle-Lomb spectrum
tauest.m:  tau estimation;
TAUrednoise.m: generate redonoise series given tau.
