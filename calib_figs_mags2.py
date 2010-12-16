import os
import numpy
import pylab
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass

filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
colors = ['b', 'g', 'y', 'r', 'm', 'k']
def color_counter_next(i):
    i = i + 1
    if i == len(colors):
        i = 0
    return i

def read_hardware(shift_perc=None):
    # read system (hardware) transmission)
    filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
    hardware = ("detector.dat", "m1.dat", "m2.dat", "m3.dat", "lens1.dat", "lens2.dat", "lens3.dat")
    # Read in the standard components, except shift the filter by 1% of the central wavelength near the edge
    sys = {}
    for f in filterlist:
        sys[f] = Bandpass()
        tlist = []
        for t in hardware:
            tlist.append(os.path.join(filterdir, t))
        sys[f].readThroughputList(tlist)
        tmpfilter = Bandpass()
        tmpfilter.readThroughput(os.path.join(filterdir, "filter_" + f + ".dat"))
        effwavelenphi, effwavelensb = tmpfilter.calcEffWavelen()
        # Add shift to wavelength of filter, if needed.
        if shift_perc != None:
            shift = effwavelensb * shift_perc
            tmpfilter.wavelen = tmpfilter.wavelen + shift
            tmpfilter.resampleBandpass()
        sys[f].wavelen, sys[f].sb = sys[f].multiplyThroughputs(tmpfilter.wavelen, tmpfilter.sb)
    return sys

def read_atmos():
    # Read some atmosphere throughputs.
    atmosdir = "."
    atmos = {}
    key = "Standard"
    atmos[key]= Bandpass()
    atmos[key].readThroughput(os.path.join(atmosdir, "atmos_std.dat"))
    key = "X=1.0" # H2O=0.7
    atmos[key] = Bandpass()
    atmos[key].readThroughput(os.path.join(atmosdir, "atmos_85233890.dat"))
    key = "X=1.5" # H2O=1.4
    atmos[key] = Bandpass()
    atmos[key].readThroughput(os.path.join(atmosdir, "atmos_85488653.dat"))
    return atmos

def combine_throughputs(atmos, sys):
    # Set up the total throughput for this system bandpass, using the variety of atmospheres. 
    total = {}
    key = 'standard'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys[f].multiplyThroughputs(atmos['Standard'].wavelen, atmos['Standard'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'X=1.0'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys[f].multiplyThroughputs(atmos['X=1.0'].wavelen, atmos['X=1.0'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'X=1.5'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys[f].multiplyThroughputs(atmos['X=1.5'].wavelen, atmos['X=1.5'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    return total


def read_seds(total):
    # read SEDs - one blue star, one red star
    stars = {}
    seddir = "/Users/rhiannonjones/seds/kurucz_r/"
    key = "red"
    stars[key] = Sed()
    stars[key].readSED_flambda(os.path.join(seddir, "km01_6000.fits_g40"))
    normmag = 20
    fluxnorm = stars[key].calcFluxNorm(normmag, total['standard']['r'])
    stars[key].multiplyFluxNorm(fluxnorm)
    stars[key].flambdaTofnu()
    key = "blue"
    stars[key] = Sed()
    stars[key].readSED_flambda(os.path.join(seddir, "kp02_35000.fits_g40"))
    normmag = 20
    fluxnorm = stars[key].calcFluxNorm(normmag, total['standard']['r'])
    stars[key].multiplyFluxNorm(fluxnorm)
    stars[key].flambdaTofnu()
    return stars



if __name__ == "__main__":
    
