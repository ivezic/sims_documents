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
    # Read in the standard components, but potentially shift the filter by shift_perc percent.
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
    for key in atmos.keys():
        total[key] = {}
        for f in filterlist:
            wavelen, sb = sys[f].multiplyThroughputs(atmos[key].wavelen, atmos[key].sb)
            total[key][f] = Bandpass(wavelen, sb)
            total[key][f].sbTophi()
    return total


def read_stars():
    # read MS, g40 stars SEDs
    stars = {}
    homedir = os.getenv("HOME")    
    stardir = os.path.join(homedir, "seds/kurucz_r")
    allfilelist = os.listdir(stardir)
    starlist = []
    for filename in allfilelist:
        if filename[-3:] == 'g40':
            starlist.append(filename)
    stars = {}
    for s in starlist:
        stars[s] = Sed()
        stars[s].readSED_flambda(os.path.join(stardir,s))
    print "#Read %d stars from %s" %(len(starlist), stardir)
    atemperature = []
    amet = []
    starlist2 = []
    for s in starlist:
        tmp = s.split('_')
        met = float(tmp[0][2:])
        if tmp[0][1] == 'm':
            met = -1 * met        
        met = met/10.0
        temperature = float(tmp[1][:5])
        if (temperature > 5000.0):
            amet.append(met)
            atemperature.append(temperature)
            starlist2.append(s)
    temperature = numpy.array(atemperature)
    met = numpy.array(amet)
    starlist = starlist2
    stars = {}
    for s in starlist:
        stars[s] = Sed()
        stars[s].readSED_flambda(os.path.join(stardir,s))
    print "#Read %d stars from %s" %(len(starlist), stardir)
    return stars, starlist, temperature, met


if __name__ == "__main__":
    
    sys_std = read_hardware(shift_perc=None)
    sys_shift = read_hardware(shift_perc=0.01)
    atmokeylist = ['standard', 'X=1.0', 'X=1.5']
    atmos = read_atmos()
    total_std = combine_throughputs(atmos, sys_std)    
    total_shift = combine_throughputs(atmos, sys_shift)
    stars, starlist, temperatures, metallicity = read_stars()
    print temperatures.max(), temperatures.min(), metallicity.max(), metallicity.min()
    print len(starlist)
    # calculate the standard and shifted magnitudes for each of these stars
    mags_std = {}
    mags_shift = {}
    for f in filterlist:
        mags_std[f] = numpy.zeros(len(starlist), dtype='float')
        mags_shift[f] = numpy.zeros(len(starlist), dtype='float')
        i = 0
        for s in starlist:
            mags_std[f][i] = stars[s].calcMag(total_std['Standard'][f])
            mags_shift[f][i] = stars[s].calcMag(total_shift['X=1.5'][f])
            i = i + 1
    gi = mags_std['g'] - mags_std['i']
    ur = mags_std['u'] - mags_std['r']
    print gi.min(), gi.max(), ur.min(), ur.max()
    shifts = {}
    for f in filterlist:
        shifts[f] = mags_shift[f] - mags_std[f]
    print len(gi), len(shifts['u'])

    pylab.figure()
    pylab.subplots_adjust(top=0.93, wspace=0.3, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
    urcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
    urbinsize = abs(ur.min() - ur.max())/6.0
    urbins = numpy.arange(ur.min(), ur.max()+urbinsize, urbinsize)
    metbinsize = abs(metallicity.min() - metallicity.max())/6.0
    metbins = numpy.arange(metallicity.min(), metallicity.max()+metbinsize, metbinsize)
    #print urbinsize, urbins
    print metbinsize, metbins
    i = 1
    for f in filterlist:
        print f
        pylab.subplot(3,2,i)
        """
        for urbidx in range(len(urbins)):
            condition = ((ur >= urbins[urbidx]) & (ur <= urbins[urbidx]+urbinsize))
            coloridx = urcolors[urbidx]
            pylab.plot(gi[condition], shifts[f][condition], coloridx+'.')
        """
        for metidx in range(len(metbins)):
            condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
            coloridx = urcolors[metidx]
            pylab.plot(gi[condition], shifts[f][condition], coloridx+'.')
        pylab.xlabel("g-i")
        pylab.ylabel("Delta %s " %(f))
        if f == 'u':
            pylab.ylim(-0.11, 0.03)
        else:
            pylab.ylim(-0.04, 0.025)
        pylab.grid(True)
        i = i + 1
    pylab.figtext(0.2, 0.95, "Change in magnitude for X=1.5 and Filter shift of 1%")
    pylab.show()
    #pylab.savefig("delta_mags2.eps", format='eps')
