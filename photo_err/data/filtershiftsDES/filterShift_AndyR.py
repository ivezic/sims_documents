"""
Read the filter shift data from AndyR and calculate a 'typical' filter shift as a function of radius.
The data from AndyR translates DES filter measurements into the LSST beam and for LSST filter edge location
specifications; there is more spatial information than simple radial terms, so more info on variation in the
focal plane could be gathered.

This is intended to be used for exploratory purposes, to understand the DESasLSSTfilters data that Andy provided.

Radius in AndyR's data is in mm in the focal plane (so, so are fits). 

"""

import os
import numpy
import pylab
import pyfits
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass

class FilterShiftDES():
    def __init__(self, dataDir=None, dataFile=None):
        """Read data from AndyR's file."""
        if dataDir == None:
            dataDir = 'data/filterShiftsDES'
        if dataFile == None:
            dataFile = 'DESfilts_as_LSSTfilts_130327.fits'
        file = os.path.join(dataDir, dataFile)
        print 'Using data file %s' %(file)
        hdus = pyfits.open(file)
        # See README in the dataDir for more info on file format.
        # Read the wavelength grid for each transmission curve (1nm resolution, 300nm-1100nm)    
        self.wavelen = hdus[1].data.field(0)[0]
        # Read the transmission curve data.
        tbdata = hdus[2].data
        # Locations of data points
        self.x = tbdata.field(0)
        self.y = tbdata.field(1)
        # r in mm
        self.rad = tbdata.field(2)
        # az in deg - r/az is same info as x/y)
        self.az = tbdata.field(3)
        # transmission function	values
        self.all = {}
        self.filterlist = ('u', 'g', 'r', 'i', 'z', 'y4')
        for i, f in enumerate(self.filterlist):
            self.all[f] = tbdata.field(i+4)
        # But note that all['u'], etc. contains 3209 points - one for each x/y or rad/az location -
        #  and at each location they then have the full transmission function.
        return

def plotR_curves(fs, az=0.0, az_range=20.0):
    """Look at throughput curves at a range of radii, to look for trends with radius. Idea being
    to understand what a 'typical' change with radius should look like. """
    # Pick curves within a particular azimuth range
    condition = (numpy.abs(fs.az - az) < az_range)
    # And then pull out the curves and tag for the radius value.
    radii = fs.rad[condition]
    wavelen = fs.wavelen
    # What we learn from this is that
    # (a) we should ignore any throughputs at higher than 500 nm in any band beyond g (artifacts)
    # (b) the peak throughput does have a trend with radius, but it's always less than 3% and generally <1.5%
    #    (b2) from looking at the results from plotF_curves we can see that these are close to, but not
    #         exactly azimuthally symmetric.
    # (c) the red and blue 50% wavelengths change differently for different filters as a function of wavelength.
    #     Generally the change (relative to the 50% wavelength at the center of the fov) is within 1% (and even
    #     less than 0.3 %). Presumably we should scale up for LSST filters / tolerances.
    #    (c2) from looking at the results from plotF_curves we can see that these are also not azimuthally symm.
    #    (c3) Fit results (using the full 360 degrees): (fit = p0*r^4 + p1*r^3 + p2*r^2 + p1*r + p0) - r in mm.
    #         (fourth order fits needed to recreate outer behavior while preserving 'flatness' near center).
    # Filter u Red Fit [  6.14143692e-11  -4.14893612e-08   9.37034194e-06  -1.12806226e-03   1.07654343e-02]
    # Filter u Blue fit [ -1.29582350e-10   8.71585854e-08  -1.66266810e-05   8.72712223e-04  -1.39733358e-02]
    # Filter g Red Fit [  5.59255744e-11  -3.08353519e-08   4.27962080e-06  -2.91183820e-04   4.22195231e-03]
    # Filter g Blue fit [  1.30441557e-10  -6.99356343e-08   8.73649283e-06  -5.62867623e-04   2.11675371e-03]
    # Filter r Red Fit [  2.61816945e-10  -1.49690259e-07   1.79901722e-05   6.26983535e-04   1.34451179e-02]
    # Filter r Blue fit [ -1.61538696e-10   3.79213416e-08   2.49058569e-05  -9.20429417e-03   8.22199892e-02]
    # Filter i Red Fit [  1.22737872e-10  -1.37037554e-07   3.54208316e-05   5.58487340e-04   5.72718292e-02]
    # Filter i Blue fit [  1.48295170e-10  -9.07946003e-08   1.05749854e-05   1.41617255e-03   4.22274038e-03]
    # Filter z Red Fit [ -3.08823985e-11  -1.13055423e-08   1.55896246e-05  -4.39261144e-03   5.10121410e-02]
    # Filter z Blue fit [ -4.53660382e-10   2.41555318e-07  -2.78940603e-05  -1.38313865e-04   6.86394071e-03]
    # Filter y4 Red Fit [  2.10879879e-11  -6.72136311e-08   3.16270858e-05  -4.44560441e-03   5.15339162e-02]
    # Filter y4 Blue fit [  4.54589095e-10  -2.45844865e-07   2.61155525e-05  -2.25361790e-05   3.88618464e-02]
    for f in fs.filterlist:
        throughputs = fs.all[f][condition]
        if f in ('r', 'i', 'z', 'y4'):
            condition2 = (wavelen < 500.)
            for i in range(len(throughputs)):
                throughputs[i][condition2] = 0.0
        # Pull out (a/the) central throughput curve
        condition2 = (radii < 0.001)
        throughputs_center = throughputs[condition2][0]
        effsb = (wavelen * throughputs_center).sum() / throughputs_center.sum()
        pylab.figure()
        pylab.title('Filter %s' %(f))
        for i in range(len(throughputs)):
            pylab.plot(wavelen, throughputs[i])
        pylab.figure()
        pylab.title('Filter %s - max throughput' %(f))
        maxthroughput= numpy.zeros(len(radii), 'float')
        for i in range(len(throughputs)):
            maxthroughput[i] = throughputs[i].max()
        pylab.plot(radii, maxthroughput, 'kx')
        pylab.xlabel('Radius')
        pylab.ylabel('Max throughput')
        red50 = numpy.zeros(len(radii), 'float')
        blue50 = numpy.zeros(len(radii), 'float')
        for i in range(len(throughputs)):
            condition2 = (wavelen < effsb)
            wavelen2 = numpy.arange(wavelen.min(), effsb, 0.02)
            tmp = numpy.interp(wavelen2, wavelen[condition2], throughputs[i][condition2])
            condition3 = (numpy.abs(tmp - throughputs[i].max()*0.5) < 0.01)
            red50[i] = wavelen2[condition3].mean()
            condition2 = (wavelen > effsb)
            wavelen2 = numpy.arange(effsb, wavelen.max(), 0.02)
            tmp = numpy.interp(wavelen2, wavelen[condition2], throughputs[i][condition2])
            condition3 = (numpy.abs(tmp - throughputs[i].max()*0.5) < 0.01)
            blue50[i] = wavelen2[condition3].mean()
        # Pull out (a/the) central red/blue edges
        condition2 = (radii < 0.001)
        red50_center = red50[condition2]
        blue50_center = blue50[condition2]
        # Calculate amount of shift (% of effsb) of red/blue edges
        rshift = (red50 - red50_center) / effsb * 100.0 
        bshift = (blue50 - blue50_center) / effsb * 100.0
        # fit red/blue edge changes
        idx = numpy.argsort(radii)
        rp = numpy.polyfit(radii[idx], rshift[idx], 4)
        ry = numpy.polyval(rp, radii[idx])
        bp = numpy.polyfit(radii[idx], bshift[idx], 4)
        by = numpy.polyval(bp, radii[idx])
        print 'Filter', f, 'Red Fit', rp
        print 'Filter', f, 'Blue fit', bp
        pylab.figure()
        pylab.title('Change in blue and red (50Perc.*peak) wavelengths %s' %(f))
        pylab.plot(radii, rshift, 'rx')
        pylab.plot(radii[idx], ry, 'k-')
        pylab.plot(radii, bshift, 'bx')
        pylab.plot(radii[idx], by, 'k-')
        pylab.xlabel('Radius')
        pylab.ylabel('Wavelength Change (Percent)')
    pylab.show()
        
def plotF_curves(fs):
    """Plot some characteristics over the 'focal plane'. """
    wavelen = fs.wavelen
    x = numpy.arange(fs.x.min(),fs.x.max())
    y = numpy.arange(fs.y.min(), fs.y.max())
    xx, yy = numpy.meshgrid(x, y)
    for f in fs.filterlist:
        throughputs = fs.all[f]
        # Max throughputs
        pylab.figure()
        pylab.title('Filter %s - max throughput' %(f))
        maxthroughput = numpy.zeros(len(fs.x), 'float')
        for i in range(len(fs.x)):
            maxthroughput[i] = throughputs[i].max()
        zz = pylab.griddata(fs.x, fs.y, maxthroughput, xx, yy, interp='linear')
        pylab.contourf(xx, yy, zz, 50)
        pylab.colorbar()
        # Red/Blue edges
        red50 = numpy.zeros(len(fs.x), 'float')
        blue50 = numpy.zeros(len(fs.x), 'float')
        for i in range(len(fs.x)):
            maxwavelen = wavelen[throughputs[i] == throughputs[i].max()]
            if len(maxwavelen) > 1:
                maxwavelen = maxwavelen.mean()
            condition2 = (wavelen < maxwavelen)
            wavelen2 = numpy.arange(wavelen.min(), maxwavelen, 0.02)
            tmp = numpy.interp(wavelen2, wavelen[condition2], throughputs[i][condition2])
            condition3 = (numpy.abs(tmp - throughputs[i].max()*0.5) < 0.01)
            red50[i] = wavelen2[condition3].mean()
            condition2 = (wavelen > maxwavelen)
            wavelen2 = numpy.arange(maxwavelen, wavelen.max(), 0.02)
            tmp = numpy.interp(wavelen2, wavelen[condition2], throughputs[i][condition2])
            condition3 = (numpy.abs(tmp - throughputs[i].max()*0.5) < 0.01)
            blue50[i] = wavelen2[condition3].mean()
        red50 = red50 / red50[0] * 100.0
        blue50 = blue50 / blue50[0] * 100.0
        pylab.figure()
        pylab.title('Change in red 50 percent wavelength %s' %(f))
        zz = pylab.griddata(fs.x, fs.y, red50, xx, yy, interp='linear')
        pylab.contourf(xx, yy, zz, 50)
        pylab.colorbar()
        pylab.figure()
        pylab.title('Change in blue 50 percent wavelength %s' %(f))
        zz = pylab.griddata(fs.x, fs.y, blue50, xx, yy, interp='linear')
        pylab.contourf(xx, yy, zz, 50)
        pylab.colorbar()
    pylab.show()
