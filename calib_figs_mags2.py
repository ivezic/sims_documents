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
    filters = read_filtersonly(shift_perc=shift_perc)
    sys = {}
    for f in filterlist:
        sys[f] = Bandpass()
        # put together the standard component list
        tlist = []
        for t in hardware:
            tlist.append(os.path.join(filterdir, t))
        # read in the standard components, combine into sys
        sys[f].readThroughputList(tlist)
        # multiply by the filter throughput for final hardware throughput (no atmosphere)
        sys[f].wavelen, sys[f].sb = sys[f].multiplyThroughputs(filters[f].wavelen, filters[f].sb)
    return sys

def read_filtersonly(shift_perc=None):
    filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
    filters = {}
    for f in filterlist:
        filters[f] = Bandpass()
        filters[f].readThroughput(os.path.join(filterdir, "filter_" + f + ".dat"))
        effwavelenphi, effwavelensb = filters[f].calcEffWavelen()
        if shift_perc != None:
            shift = effwavelensb * shift_perc/100.0
            print f, shift
            filters[f].wavelen = filters[f].wavelen + shift
            filters[f].resampleBandpass()
    return filters
        
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

    # shift the filters by nothing (standard)
    sys_std = read_hardware(shift_perc=None)
    # shift the filters by one percent
    shift_perc = 0.05
    sys_shift = read_hardware(shift_perc=shift_perc)
    # read in a few different atmospheres
    atmokeylist = ['standard', 'X=1.0', 'X=1.5']
    atmos = read_atmos()
    # combine to total throughputs 
    total_std = combine_throughputs(atmos, sys_std)    
    total_shift = combine_throughputs(atmos, sys_shift)
    # read the stars
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
            #mags_shift[f][i] = stars[s].calcMag(total_shift['X=1.5'][f])
            atmo_choice = 'Standard'  # atmo='X=1.5'
            mags_shift[f][i] = stars[s].calcMag(total_shift[atmo_choice][f])
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
            pylab.ylim(-0.005, 0.0025)
        else:
            #pylab.ylim(-0.04, 0.025)
            pylab.ylim(-0.0025, 0.0025)
        pylab.grid(True)
        i = i + 1
    pylab.figtext(0.2, 0.95, "Change in magnitude for %s atmosphere and Filter shift of %.2f%s" %(atmo_choice, shift_perc, "%"))
    #pylab.savefig("delta_mags2.eps", format='eps')

    if True:
        # plot the shift in the filters 
        pylab.figure()
        filters_std = read_filtersonly()
        filters_shift = read_filtersonly(shift_perc=shift_perc)
        # plot the filters alone
        pylab.subplot(211)
        i = 0
        for f in filterlist:
            pylab.plot(filters_std[f].wavelen, filters_std[f].sb, colors[i]+"-", label=f)
            pylab.plot(filters_shift[f].wavelen, filters_shift[f].sb, colors[i]+":")
            i = color_counter_next(i)
        pylab.ylim(0, 1)
        pylab.xlim(300, 1200)
        pylab.ylabel("Transmission")
        pylab.grid(True)
        # plot the total throughput
        pylab.subplot(212)
        i = 0
        for f in filterlist:
            pylab.plot(total_std['Standard'][f].wavelen, total_std['Standard'][f].sb, colors[i]+"-", label=f)
            pylab.plot(total_shift[atmo_choice][f].wavelen, total_shift[atmo_choice][f].sb, colors[i]+":")
            i = color_counter_next(i)
        pylab.ylim(0, 0.8)
        pylab.xlim(300, 1200)
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Transmission")
        pylab.grid(True)
        pylab.legend(numpoints=1, fancybox=True, loc='upper right')
            
    pylab.show()

