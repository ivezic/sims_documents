"""
FilterShifts:
a class to implement filter shifts as fit from AndyR's data (where he 'interpreted' DES filters as LSST
filters - by altering the throughput curves for our different incident angle range as well as different 50%
points and bandpass edge requirements).
The filter shifts from the interpreted DES filters are generally small - <1% * effsb - so we also (optionally)
scale these overall shifts up to the levels permitted by the LSST requirements (<2.5% * effsb).

These shifts can then be used to calculate the potential errors in determining color-dependent magnitude
corrections due to the bandpass shift --- these errors could come from (a) assuming the bandpass was constant
over the entire field of view [although this should not be the case, since we will have the monochromatic
flat fields to compensate for this .. the monochromatic flats will have errors which can be predicted using the
GhostData class] or (b) from filter 'jitter' where the filter does not return to the exact same position it had
during the monochromatic flat measurements [these errors are unavoidable].

The color-dependent magnitude errors are calculated assuming that we have errors due to (b), rather than (a) ..
although you could easily calculate what the overall color-terms due to the radial change in bandpass is, as well.

### Waiting for sysEng to confirm max_jitter requirements and proper max scale values. 

"""

import os
import numpy
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass


# Red and Blue shifts fit from DES filters interpreted to LSST specs (as per AndyR and fit using the code
#  in data/filterShift_AndyR.py). These values are numpy.polyfit from
#        rshift = (red50 - red50_center) / effsb * 100.0 
#        bshift = (blue50 - blue50_center) / effsb * 100.0
#  (fit = p0*r^4 + p1*r^3 + p2*r^2 + p1*r + p0) - r in mm.
# Scale up total shift from DES-interpreted values to values allowed by requirements (i.e. from <1% to 2.5% max)
#   using shift_scale (use same scale for both, because otherwise end up pushing both 'edges' to the same change,
#   which I don't think we want .. want to simulate 'stretch' of bandpass due to difference in changes in blue/red
#   edges, if these differences exist).  Note that you could potentially re-scale 'shift_scale' by remembering
#   that it was calculated to make the max shift = 2.5% ... so if you wanted max-shift=1%, just /2.5, etc.
red_shift = {'u': [6.14143692e-11, -4.14893612e-08, 9.37034194e-06, -1.12806226e-03, 1.07654343e-02],
             'g': [5.59255744e-11, -3.08353519e-08,  4.27962080e-06, -2.91183820e-04,  4.22195231e-03],
             'r': [2.61816945e-10, -1.49690259e-07,  1.79901722e-05,  6.26983535e-04,  1.34451179e-02],
             'i': [1.22737872e-10, -1.37037554e-07,  3.54208316e-05,  5.58487340e-04,  5.72718292e-02],
             'z': [-3.08823985e-11, -1.13055423e-08,  1.55896246e-05, -4.39261144e-03,  5.10121410e-02],
             'y4': [2.10879879e-11, -6.72136311e-08,  3.16270858e-05, -4.44560441e-03,  5.15339162e-02],
             'y': [2.10879879e-11, -6.72136311e-08,  3.16270858e-05, -4.44560441e-03,  5.15339162e-02]}
blue_shift = {'u':[-1.29582350e-10,  8.71585854e-08, -1.66266810e-05,  8.72712223e-04, -1.39733358e-02],
              'g':[1.30441557e-10, -6.99356343e-08,  8.73649283e-06, -5.62867623e-04,  2.11675371e-03],
              'r':[ -1.61538696e-10,  3.79213416e-08,  2.49058569e-05, -9.20429417e-03,  8.22199892e-02],
              'i':[1.48295170e-10, -9.07946003e-08, 1.05749854e-05, 1.41617255e-03,  4.22274038e-03],
              'z':[-4.53660382e-10,  2.41555318e-07, -2.78940603e-05, -1.38313865e-04,  6.86394071e-03],
              'y4':[4.54589095e-10, -2.45844865e-07,  2.61155525e-05, -2.25361790e-05,  3.88618464e-02],
              'y':[4.54589095e-10, -2.45844865e-07,  2.61155525e-05, -2.25361790e-05,  3.88618464e-02]}
shift_scale = {'u': 21.167527337, 'g': 11.6569392979, 'r': 2.58212509247, 'i': 3.61083558372, 
               'z': 4.34369995266, 'y4': 3.6532190412, 'y': 3.6532190412 }
# Radius range for the original fits (AndyR's data)
RADIUS_MIN = 0.0   #mm
RADIUS_MAX = 350.0 #mm

class FilterShift():
    def __init__(self, throughputsDir=None, filterlist=('u', 'g', 'r', 'i', 'z', 'y4')):
        """Instantiate object and do a reasonable default setup ready for calculating magnitudes."""
        self.read_base_throughputs(throughputsDir, filterlist)
        self.setBandpasses(max_jitter=2.0)
        self.setPhiArray()
        return

    def read_base_throughputs(self, throughputsDir=None, filterlist=('u', 'g', 'r', 'i', 'z', 'y4')):
        """Read base throughput curves from throughputsDir (if None, uses throughputs package default).
        Saves filter throughput curves separately from atmosphere/mirrors/lenses/detector, to enable shifting."""
        self.filterlist = filterlist
        # Read the basic default throughput curves (except for filters).
        if throughputsDir == None:
            throughputsDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')        
        components = ('atmos.dat', 'm1.dat', 'm2.dat', 'm3.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat',
                      'detector.dat')
        self.base = Bandpass()
        self.base.readThroughputList(componentList=components, rootDir=throughputsDir)
        self.base_filters = {}
        for f in self.filterlist:
            self.base_filters[f] = Bandpass()
            # Read the basic filter throughput curve. 
            self.base_filters[f].readThroughput(os.path.join(throughputsDir, 'filter_' + f + '.dat'))
        print '# Read base throughput curves from %s' %(throughputsDir)
        return

    def combine_throughputs(self, radius=0):
        """Combine base throughput curves from atmosphere/mirrors/lenses/detector with shifted filter
        throughput curves at a particular radius. Calls shift_filter to shift the filter. """
        # Combine throughputs from base and from the filters (at a particular radius) to make a full
        #  throughput curve.
        # Check to see if overall (complete) bandpass dictionary exists.
        try:
            self.bp
        except:
            self.bp = {}
        # Set up a dictionary for this radius, unless it already exists (in which case, we're done).
        try:
            self.bp[radius][self.filterlist[0]]
            return
        except:
            self.bp[radius] = {}
        # Calculate combined throughput curves at this radius. 
        for f in self.filterlist:
            self.bp[radius][f] = Bandpass()
            # Shift filter to appropriate values at this radius.
            fwavelen, fsb = self.shift_filter(f, radius)
            # Note that multiplyThroughputs will resample fwavelen/fsb onto the same grid as base. 
            wavelen, sb = self.base.multiplyThroughputs(fwavelen, fsb)
            # Set bp[r][f].
            self.bp[radius][f].setBandpass(wavelen, sb)
            self.bp[radius][f].sbTophi()
        return

    def shift_filter(self, f, radius, scale=True):
        """Shift/stretch the filter red and blue edges by the desired amount, for this filter f
        and at radius 'radius', and return wavelen/sb. """
        # Calculate red and blue shifts at this radius (radius must be a single value).
        #   numpy.polyval(r_shift[f]) gives the shift % = (red50 - red50_center) / effsb * 100.0 
        # and then this translates into an actual value to add to the red wavelengths as
        #    (%/100.*effsb) = red50 - red50_baseline.   (red50 = red50_baseline + shift/100.*effsb)
        # This will also be scaled up to LSST permitted shift values, if scale=True. (otherwise max shift <.5%). 
        rshift = numpy.polyval(red_shift[f], radius)
        bshift = numpy.polyval(blue_shift[f], radius)
        if scale==True:
            rshift = rshift * shift_scale[f]
            bshift = bshift * shift_scale[f]
        # Because we have different shifts on blue/red edges, split at effsb and stretch each side.
        effsb = self.base_filters[f].calcEffWavelen()[1]
        wavelen = numpy.copy(self.base_filters[f].wavelen)
        # Shift the red side
        condition = (wavelen > effsb)
        wavelen[condition] = wavelen[condition] + rshift / 100.0 * effsb
        # Shift the blue side
        condition = (wavelen < effsb)
        wavelen[condition] = wavelen[condition] + bshift / 100.0 * effsb
        # Wavelen now represents the shifted bandpass (using the original throughput values, but 'stretched'). 
        return wavelen, self.base_filters[f].sb
    
    def setBandpasses(self, max_jitter=2.0, radius_min=RADIUS_MIN, radius_max=RADIUS_MAX):
        """Set up Bandpass objects covering the radius range with steps of max_jitter, for each filter."""
        # We must compare mags for shifted (at radius 'r') bandpass and mags at the same radius but for a filter with a 
        #  'jitter' in its position. The max jitter (assume = max error) is equivalent to looking at a radius +/- the max jitter amount.
        # Set these up for a series of radii, separated by max jitter amount.
        self.radii = numpy.arange(radius_min, radius_max+max_jitter, max_jitter)
        for r in self.radii:
            # Generate self.bp[r][f]
            self.combine_throughputs(r)
        return

    def setPhiArray(self):
        """After setting up Bandpass objects over radius range, generate phi arrays for use in manyMagCalc. """
        # Generate phiArrays  (one per filter, covering the whole radius range).
        self.phiarray = {}
        for f in self.filterlist:
            self.phiarray[f] = numpy.empty((len(self.radii), len(self.base.wavelen)), dtype='float')
            # Go through radii and generate combined (with shifted filter) throughputs, add to phiArray.
            for i, r in enumerate(self.radii):
                self.phiarray[f][i] = self.bp[r][f].phi
        self.wavelen_step = self.base.wavelen[1] - self.base.wavelen[0]
        return
    
    def plotBandpasses(self):
        """Just to provide a sanity check and visual inspection of some of the shifted bandpasses. """
        import pylab
        r0 = self.radii[0]
        r1 = self.radii[len(self.radii)/2]
        r2 = self.radii[len(self.radii)-1]
        for f in self.filterlist:
            pylab.figure()
            pylab.plot(self.bp[r0][f].wavelen, self.bp[r0][f].sb, label='Center')
            pylab.plot(self.bp[r1][f].wavelen, self.bp[r1][f].sb, label='Middle')
            pylab.plot(self.bp[r2][f].wavelen, self.bp[r2][f].sb, label='Edge')
            pylab.xlabel('Wavelength (nm)')
            pylab.ylabel('Throughput (0-1)')
            pylab.title('Filter Shift for %s' %(f))
            pylab.legend(fontsize='smaller', fancybox=True)
        pylab.show()
        return
    
    def calc_mags(self, sed, f):
        """Calculate magnitudes for 'sed' in filter 'f', at over the radii set up in setPhiArray. 
        Returns mags at each radii, as well as dmags representing the max error due to filter jitter (=steps in radius array), 
         (equivalent to the max difference between +/- one step in mag). """
        # Make sure that sed's fnu exists and is on the same wavelength grid as the phiarray grid.
        r0 = self.radii[0]
        wavelen_match = self.base.wavelen
        if sed.needResample(wavelen_match=wavelen_match):
            sed.resampleSED(wavelen_match=wavelen_match)
        sed.flambdaTofnu()
        # Calculate the magnitudes for the bandpass as would be measured (i.e. @ radius, we're not including ghosting induced errors) 
        #  and as might be the result with jitter. Assuming max error happens in the max jitter 'direction', this means looking at radii
        #  at values +/- jitter to look for the max difference in magnitude. 
        mags = sed.manyMagCalc(self.phiarray[f], self.wavelen_step)
        # And the color-dependent differences in natural magnitudes.
        dmags_up = mags[:-1] - mags[1:]
        dmags_up = numpy.concatenate((dmags_up, [0]))
        dmags_down = mags[1:] - mags[:-1]
        dmags_down = numpy.concatenate(([0], dmags_down))
        # Return the value of the largest offset (absolute value) at each radius, in mmags.
        dmags = numpy.where(numpy.abs(dmags_up)>numpy.abs(dmags_down), dmags_up, dmags_down)
        dmags = dmags * 1000.0
        # yes, those steps above are a bit overkill/unnecessary ... but it does keep absolutely straight the radius/dmags relationship. 
        return mags, dmags
    

