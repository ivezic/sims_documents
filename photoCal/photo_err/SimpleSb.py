"""
SimpleSb - a class to allow simple changes to Sb (and thus phi) to evaluate effects of
  - shifting bandpass (by given amount, at all wavelenghts)
  - changing resolution of measuring bandpass
  - ?? (TBD)
to help evaluate requirements on phi for stars/SN/emission line galaxies. 
Calculate mag(orig) / mag(changed) / dmag (=orig-changed)  (= color-term errors). 

Comment: final resolution of Bandpass must be < 1 Angstrom in order to be sure to resolve emission lines in emission line galaxy. (since Sed/Bandpass
share the same resolution when actually calculating the magnitudes). 

Class uses basic (baseline) total throughput curves to operate on. (or substitute your own curves in readBandpasses). 
"""

import os
from copy import deepcopy
import numpy
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass


class SimpleSb():
    def __init__(self, wavelen_step=0.05, throughputsDir=None, filterlist=('u', 'g', 'r', 'i', 'z', 'y'), rootName='total_', suffixName='.dat'):
        """Instantiate the object and read base throughput curves. """
        self.wavelen_step = wavelen_step # in nm
        self.read_throughputs(throughputsDir=throughputsDir, filterlist=filterlist, rootName=rootName, suffixName=suffixName)
        return

    def read_throughputs(self, throughputsDir=None, filterlist=('u', 'g' , 'r', 'i', 'z', 'y'), rootName='total_', suffixName='.dat'):
        """Read total throughputs curves for filterlist from throughputsDir, with rootnames 'rootName' (suffix='suffixName').  """
        # Set throughputsDir if not given.
        if throughputsDir == None:
            throughputsDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
        self.filterlist = filterlist
        # Read in total throughputs curves, save in dictionary 'basebp'.
        self.basebp = {}
        for f in filterlist:
            self.basebp[f] = Bandpass()
            self.basebp[f].readThroughput(os.path.join(throughputsDir, rootName + f + suffixName), wavelen_step=self.wavelen_step)
            self.basebp[f].sbTophi()
        return

    def shift_throughputs(self, shift):
        """Shift base throughputs curves by amount shift. 
        If 'shift' is a dictionary (can specify different shifts for different filters, but keys must match self.filterlist), 
        then will be applied as such (otherwise same shift applied to all filters). """
        try:
            self.newbp 
            print '# self.newbp already exists - using existing newbp and adding this shift.'
        except:
            self.newbp = deepcopy(self.basebp)
            print '# creating newbp and applying shift.'
        for f in self.filterlist:
            # Shift each bandpass.
            if isinstance(shift, dict):
                self.newbp[f].wavelen += shift[f]
            else:
                self.newbp[f].wavelen += shift
            # Resample onto the original grid (so that wavelengths start/end at intended locations).
            self.newbp[f].resampleBandpass()
            self.newbp[f].sbTophi()
        return


    def deres_throughputs(self, resolution):
        """Sample bandpass with steps corresponding to 'resolution' (which can be a dictionary, per filter, or single #), 
        then do linear interpolation between these sample points to regenerate a bandpass with original wavelength spacing. """
        try:
            self.newbp
            print '# self.newbp already exists - using existing newbp and degrading resolution.'
        except:
            self.newbp = deepcopy(self.basebp)
            print '# creating newbp and degrading resolution.'
        for f in self.filterlist:
            if isinstance(resolution, dict):
                res = resolution[f]
            else:
                res = resolution
            # Sample bandpass with resolution 'res'
            wavelen = numpy.arange(self.newbp[f].wavelen_min, self.newbp[f].wavelen_max+res, res)
            sb = numpy.interp(wavelen, self.newbp[f].wavelen, self.newbp[f].sb, left=0., right=0.)
            # And then generate a new (linear interpolation) of bandpass with original wavelength spacing.
            self.newbp[f].setBandpass(wavelen, sb, wavelen_step = self.wavelen_step)
            self.newbp[f].sbTophi()
        return

    def plot_throughputs(self):
        """Sanity check and visualize bandpasses - original (base) and new (newbp). """
        import pylab
        for f in self.filterlist:
            pylab.figure()
            pylab.plot(self.basebp[f].wavelen, self.basebp[f].sb, label='Base')
            pylab.plot(self.newbp[f].wavelen, self.newbp[f].sb, label='New')
            pylab.legend(fancybox=True, numpoints=1)
            pylab.xlabel('Wavelength (nm)')
            pylab.ylabel('Throughput (0-1)')
            pylab.title('Filter %s' %(f))
        return

    def plot_throughputs_diff(self):
        """Sanity check and visualize bandpasses - original (base) and new (newbp). """
        import pylab
        pylab.figure()
        for f in self.filterlist:
            residuals = self.basebp[f].phi-self.newbp[f].phi
            idx = numpy.arange(numpy.size(residuals))#numpy.where(numpy.abs(residuals) > .001)
            pylab.plot(self.basebp[f].wavelen[idx], \
                       residuals[idx], label=f)
            pylab.xlabel('Wavelength (nm)')
            pylab.ylabel(r'$\Delta$ $\Phi$ (nm$^{-1}$)')
            pylab.legend(loc=4)
        return
    
    def calc_mags(self, sed, f):
        """Calculate (natural) magnitudes in original and new bandpasses, to evaluate colorterm errors induced by changes in phi. 
        Returns mag in original bandpass, mag in new bandpass (in magnitudes), and dmag (orig-new) in mmag. """
        mag_orig = sed.calcMag(self.basebp[f])
        mag_new = sed.calcMag(self.newbp[f])
        dmag = (mag_orig - mag_new) * 1000.0
        return mag_orig, mag_new, dmag


        
