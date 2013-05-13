"""
GhostData : a class to read AndyR's ghost data files (and other throughput files that contribute to the total
 throughput of the LSST system) and hold this data in
   self.direct[filter][radius]
   self.flatghost[filter][radius]
where radius == the radii where AndyR calculated the direct (primary light path only),
and flatghost (the direct + ghosting light == the light that would be observed in a uniformly illuminated flat
field image, could also be called 'flatfield'), and filter = the various LSST filters.

This can then be used to estimate the errors produced by our imperfect knowlege of the bandpass due to ghosting. 

After checking with AndyR, I found his data includes the rest of the *camera* hardware system
(detector + lenses) but not the atmosphere or mirrors.
Which version, exactly, of the detector + lenses was used is not entirely clear yet, but
will add more info when available.
This means that calculating the error due to misrepresenting the bandpass, we still need to add in the atmosphere
and mirrors when calculating the dmag(color_term) for stars, but do not need to include the other
throughput components.
When calculating the gray-scale IC for the flatfield only, then should not include the atmosphere
(but should include the mirrors), for the flat-field measurement (i.e. in flatghost) but should include atmosphere
as well as mirrors in the 'direct' path, as stars (real photometric measurements) will include the atmosphere.

Requires pyfits and the Bandpass class (from catalogs_measures), as well as an LSST throughputs set, such as from
the throughputs package.

Calculates the color-dependent magnitude changes for a given SED,
 as well as the gray-scale illumination correction.

Dmags =  true - wrong.
"""

import os 
import numpy
import pyfits
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass

class GhostData():
    def __init__(self, ghostData=None, ghostDataDir=None, throughputsDir=None,
                 filterlist=('u', 'g', 'r', 'i', 'z', 'y4')):
        """Instantiate object and do setup for magnitude calculations. """
        self.read_ghosting(ghostData, ghostDataDir, throughputsDir, filterlist)
        return

    def read_ghosting(self, ghostData=None, ghostDataDir=None,
                      throughputsDir=None, filterlist=('u', 'g', 'r', 'i', 'z', 'y4')):
        """Read in the ghost file from AndyR's fits file, and read the throughputs files to combine with
        his (filter-only) data."""
        # Read in AndyR's data from the fits file he provides (for various vendors).
        # Open and read the (various vendor) ghost data files produced by AndyR.
        if ghostData == None:
            ghostData = 'camera_ghosting_ff_calibbias_jdsu_lsstcone_121128.fits'
        if ghostDataDir == None:
            ghostDataDir = 'data/ghostingAndyR'
        file = os.path.join(ghostDataDir, ghostData)
        print '# Using data file from %s' %(file)
        hdus = pyfits.open(file)
        # Andy has always put these filters in extensions (1, 2, 3, 4, 5, 6) respectively,
        #  although in general this should be checked with a match against hdus.info()
        direct_only = {}
        flatghost_only = {}
        # Read information from file, reorganize. 
        for i, f in enumerate(filterlist):
            direct_only[f] = {}
            flatghost_only[f]= {}
            tbdata = hdus[i+1].data
            # Then order of this data is radius, wavelength, direct, flatghost [0, 1, 2, 3]
            radii = numpy.unique(tbdata.field(0))
            condition = (tbdata.field(0) == radii[0])
            wavelen = tbdata.field(1)[condition]
            for r in radii:
                # Read direct throughput data from AndyR file
                condition = (tbdata.field(0) == r)
                tmp = tbdata.field(2)[condition]
                direct_only[f][r] = Bandpass()
                direct_only[f][r].setBandpass(wavelen=wavelen, sb=tmp)
                # Resample so that we know it's on 'standard' wavelength grid we use for throughputs.
                direct_only[f][r].resampleBandpass()
                direct_only[f][r].sbTophi()
                # Read ghost+direct throughput data from AndyR file
                tmp = tbdata.field(3)[condition]
                flatghost_only[f][r] = Bandpass()
                flatghost_only[f][r].setBandpass(wavelen=wavelen, sb=tmp)
                flatghost_only[f][r].resampleBandpass()
                flatghost_only[f][r].sbTophi()
            del tbdata
        hdus.close()
        # Now build throughput curves we need for calculating gray scale and color-terms.
        # Add mirror throughputs to everything, add atmosphere throughputs to color-terms and
        #   direct throughputs used in gray IC calculation (i.e. direct has only one version, flatghost has 2)
        # Read in the basic hardware and atmosphere
        if throughputsDir == None:
            throughputsDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
        #components = ('atmos.dat', 'm1.dat', 'm2.dat', 'm3.dat',
        #              'lens1.dat', 'lens2.dat', 'lens3.dat', 'detector.dat')
        components = ('atmos.dat', 'm1.dat', 'm2.dat', 'm3.dat')
        base1 = Bandpass()
        base1.readThroughputList(componentList=components, rootDir=throughputsDir)
        components = ('m1.dat', 'm2.dat', 'm3.dat')
        base2 = Bandpass()
        base2.readThroughputList(componentList=components, rootDir=throughputsDir)
        direct = {}
        flatghost = {}
        flatghost_gray = {}
        for f in filterlist:
            direct[f] = {}
            flatghost[f] = {}
            flatghost_gray[f] = {}
            for r in radii:
                # Combine base1 with AndyR's ghosting data for color-terms.
                wavelen, tmp = base1.multiplyThroughputs(direct_only[f][r].wavelen, direct_only[f][r].sb)
                direct[f][r] = Bandpass()
                direct[f][r].setBandpass(wavelen=wavelen, sb=tmp)
                direct[f][r].sbTophi()
                wavelen, tmp = base1.multiplyThroughputs(flatghost_only[f][r].wavelen, flatghost_only[f][r].sb)
                flatghost[f][r] = Bandpass()
                flatghost[f][r].setBandpass(wavelen=wavelen, sb=tmp)
                flatghost[f][r].sbTophi()
                # Combine base2 with AndyR's ghosting data for flatghost (direct+ghost) gray-scale IC calculation.
                wavelen, tmp = base2.multiplyThroughputs(flatghost_only[f][r].wavelen, flatghost_only[f][r].sb)
                flatghost_gray[f][r] = Bandpass()
                flatghost_gray[f][r].setBandpass(wavelen=wavelen, sb=tmp)
                flatghost_gray[f][r].sbTophi()
        self.direct_only = direct_only
        self.flatghost_only = flatghost_only
        self.direct = direct
        self.flatghost = flatghost
        self.flatghost_gray = flatghost_gray
        self.radii = radii
        self.filterlist = filterlist
        self.setupPhiArray()
        return

    def setupPhiArray(self):
        """Set up phi arrays for fast magnitude calculations at all radii. """
        # Set up phi arrays (for fast mag calculations).
        # Set up phi arrays
        self.direct_phiarray = {}
        self.flatghost_phiarray = {}
        self.wavelen_step  = {}
        r0 = self.radii[0]
        for f in self.filterlist:
            self.direct_phiarray[f]=numpy.empty((len(self.radii), len(self.direct[f][r0].wavelen)), dtype='float')
            self.flatghost_phiarray[f]=numpy.empty((len(self.radii), len(self.direct[f][r0].wavelen)),
                                                   dtype='float')
            self.wavelen_step[f] = self.direct[f][r0].wavelen[1] - self.direct[f][r0].wavelen[0]
            for i, r in enumerate(self.radii):
                self.direct_phiarray[f][i] = self.direct[f][r].phi
                self.flatghost_phiarray[f][i] = self.flatghost[f][r].phi
        return

    def apply_CCOB(self, percent_reduction=90):
        """Apply a correction to the flatghost to remove percent_reduction of the ghost (bringing it closer to 'true' / direct). 
        Percent_reduction is applied to remove (percent_reduction) of the difference between flatghost - direct. """
        # This should perhaps be applied to difference of flatghost_only/direct_only, but going with entire term for now. 
        # (or maybe difference of these? need clarification from Chuck)
        for f in self.filterlist:
            for r in self.radii:
                diff = self.flatghost[f][r] - self.direct[f][r]
                reduced_diff = diff * percent_reduction / 100.0
                self.flatghost[f][r] = self.direct[f][r] + reduced_diff
        self.setupPhiArray()
        return                        

    def plot_ghosting(self, f, vendor='', xlim=[300, 1100]):
        """Plot the direct vs. direct+ghosting wavelength response."""
        import pylab
        r0 = self.radii[0]
        rn = self.radii[len(self.radii)-1]
        pylab.figure()
        pylab.subplot(211)
        pylab.title('Direct(black) and Flatghost(red) throughputs, %s %s band\n Center' %(vendor, f))
        pylab.plot(self.direct[f][r0].wavelen, self.direct[f][r0].sb, 'k-')
        pylab.plot(self.flatghost[f][r0].wavelen, self.flatghost[f][r0].sb, 'r-')
        pylab.xlim(xlim[0], xlim[1])
        pylab.ylabel('Throughput, 0-1')
        pylab.subplot(212)
        eps = 1e-21
        ratio = numpy.where(self.direct[f][r0].sb>eps, self.flatghost[f][r0].sb / self.direct[f][r0].sb, 1)
        pylab.plot(self.direct[f][r0].wavelen, ratio, 'b-')
        pylab.xlim(xlim[0], xlim[1])
        pylab.xlabel('Wavelength (nm)')
        pylab.ylabel('Flatghost/Direct')
        pylab.figure()
        pylab.subplot(211)
        pylab.title('Direct(black) and Flatghost(red) throughputs, %s %s band\n Edge' %(vendor, f))
        pylab.plot(self.direct[f][rn].wavelen, self.direct[f][rn].sb, 'k-')
        pylab.plot(self.flatghost[f][rn].wavelen, self.flatghost[f][rn].sb, 'r-')
        pylab.xlim(xlim[0], xlim[1])
        pylab.ylabel('Throughput, 0-1')
        pylab.subplot(212)
        ratio = numpy.where(self.direct[f][rn].sb>eps, self.flatghost[f][rn].sb / self.direct[f][rn].sb, 1)
        pylab.plot(self.direct[f][rn].wavelen, ratio, 'b-')
        pylab.xlim(xlim[0], xlim[1])
        pylab.xlabel('Wavelength (nm)')
        pylab.ylabel('Flatghost/Direct')
        pylab.figure()
        pylab.title('Difference between Direct and Flatghost Phi, %s %s band\n Center (Blue) vs Edge (Red)'
                    %(vendor, f))
        pylab.plot(self.direct[f][rn].wavelen, self.direct[f][r0].phi - self.flatghost[f][r0].phi, 'b-')
        pylab.plot(self.flatghost[f][rn].wavelen, self.direct[f][rn].phi - self.flatghost[f][rn].phi, 'r-')
        pylab.xlabel('Wavelength (nm)')
        pylab.ylabel(r'$\Delta Phi')
        pylab.show()
        return

    def calc_mags(self, sed, f):
        """Calculate magnitudes (note: mags with AB zeropoint & calculated using phi(lambda) -- thus, delta(mag)
        means a color-dependent delta mag only, not gray-scale!) for a particular SED, in both direct and
        direct+ghosting (flatghost) bandpasses. In filter 'f', and at all radii.
        Returns direct magnitudes (no ghosting), flatghost mags (natural mag we would assume the object had, given
        the bandpass as measured with ghosting by the monochromatic flat), and dmags (flatghost-direct, in mmags).
        Note that order of these arrays is the same as self.radii."""
        # Make sure that sed's fnu exists and is on the same wavelength grid as the phiarray grid.
        r0 = self.radii[0]
        wavelen_match = self.direct[f][r0].wavelen
        if sed.needResample(wavelen_match=wavelen_match):
            sed.resampleSED(wavelen_match=wavelen_match)
        sed.flambdaTofnu()
        # Calculate the magnitudes. 
        direct_mags = sed.manyMagCalc(self.direct_phiarray[f], self.wavelen_step[f])
        flatghost_mags = sed.manyMagCalc(self.flatghost_phiarray[f], self.wavelen_step[f])
        # And the color-dependent differences in natural magnitudes.
        dmags = (direct_mags - flatghost_mags) * 1000.0
        return direct_mags, flatghost_mags, dmags
    
    def calc_grayscale(self, f):
        """Calculate the gray-scale change for a flat SED, for filter 'f' and at all radii.
        This is essentially calculating the gray-scale illumination correction.
        Returns dmags (flatghost - direct, in mmags). """
        ic_gray = numpy.zeros(len(self.radii), 'float')
        for i, rad in enumerate(self.radii):
            ic_gray[i] = ((-2.5*numpy.log10(self.direct[f][rad].sb.sum()*self.wavelen_step[f])) - 
                          (-2.5*numpy.log10(self.flatghost_gray[f][rad].sb.sum()*self.wavelen_step[f])))
        ic_gray = (ic_gray - ic_gray.min()) * 1000.0
        return ic_gray
