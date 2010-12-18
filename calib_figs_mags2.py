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


def read_stars():
    # read MS, g40 stars SEDs
    stars = {}
    stardir = "/Users/rhiannonjones/seds/kurucz_r/"
    allfilelist = os.listdir(stardir)
    starlist = []
    for filename in allfilelist:
        if filename[-3:] == 'g40':
            starlist.append(filename)
    stars = {}
    for s in starlist:
        stars[s] = Sed.Sed()
        stars[s].readSED_flambda(stardir+s)
    print "#Read %d stars from %s" %(len(starlist), stardir)
    atemperature = []
    amet = []
    for s in starlist:
        tmp = s.split('_')
        met = float(tmp[0][2:])
        if tmp[0][1] == 'm':
            met = -1 * met
        met = met/10.0
        amet.append(met)
        temperature = float(tmp[1][:5])
        atemperature.append(temperature)
    temperature = n.array(atemperature)
    met = n.array(amet)
    return stars, starlist, temperature, met


if __name__ == "__main__":
    
    sys_std = read_hardware(shift_perc=None)
    atmos = read_atmos()
    total_std = combine_throughputs(atmos, sys_std)    
    stars, starlist = read_stars()
    # Use a bandpass to set up the fnuArray for easy/fast calculation
    # of magnitudes for the stars.
    starfnuArray, starzp = total['standard']['g'].setupFnuArray(stars)
    # calculate the standard magnitudes for each of these stars
    mags_std = {}
    atmokeylist = ['standard', 'X=1.0', 'X=1.5']
    for akey in atmokeylist:
        bplist = []
        for f in filterlist:
            bplist.append(total[akey][f])

    
    shifts = numpy.arange(0, 10, 0.5)
    for shift in shifts:
        sys = read_hardware(shift_perc = shift)
        total = combine_throughputs(atmos, sys)
        
        # calculate magnitude in flat bandpass
        for s in starlist:
            stdmags = {}
        writestring = "%s  standard " %(s)
        for f in filterlist:
            stdmags[f] = stars[s].calcMag(standard[f])
            writestring = writestring + " %.5f " %(stdmags[f])
        for f in filterlist:
            writestring = writestring + " 0.0 "
        print writestring
        for a in atmoslist:
            writestring = "%s  %s " %(s, a)
            mags = {}
            for f in filterlist:
                mags[f] = stars[s].calcMag(throughputs[a][f])
                writestring = writestring + " %.5f " %(mags[f])
            for f in filterlist:
                tmp = mags[f] - stdmags[f]
                writestring = writestring + " %.5f " %(tmp)
            print writestring

