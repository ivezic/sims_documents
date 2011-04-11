import os
import numpy
import pylab
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import AtmoComp as ac

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

###
def get_atmosDict():
    atmos = {}
    atmos['Standard'] = read_stdatmo()
    atmocmp = ac.AtmoComp()
    # max atmo = update_atmos(atmocmp, X=1.2, t0=5.6/100.0, alpha=-1.8, O3=1.5, H2O=1.3)
    # min atmo = update_atmos(atmocmp, X=1.2, t0=0.2/100.0, alpha=-0.5, O3=0.6, H2O=0.5)
    # 30p atmo = update_atmos(atmocmp, X=2.5, t0=(0.8/100), alpha=-1.0, O3=0.9, H2O=0.8)
    # 30p atmo = update_atmos(atmocmp, X=2.5, t0=(2.4/100.0), alpha=-1.4, O3=1.17, H2O=1.04)
    # 10p/30 atmo = update_atmos(atmocmp, X=1.2, t0=(0.8/100.0), alpha=-1.0, O3=0.9, H2O=0.8)
    # 10p/30 atmo = update_atmos(atmocmp, X=1.2, t0=(1.3/100.0), alpha=-1.13, O3=0.99, H2O=1.04)
    for X in ('1.2', '2.5'):
        #update_atmos(atmocmp, X=float(X))
        if X=='1.2':
            atmo = update_atmos(atmocmp, X=2.0, t0=(0.8/100.0), alpha=-1.0, O3=0.9, H2O=0.8)
        if X=='2.5':
            atmo = update_atmos(atmocmp, X=2.0, t0=(1.3/100.0), alpha=-1.13, O3=0.99, H2O=1.04)
        atmos[X] = atmo_BP(atmocmp)
    return atmos

def read_stdatmo():
    atmosdir = "."
    atmos_bp = Bandpass()
    atmos_bp.readThroughput(os.path.join(atmosdir, "atmos_std.dat"))
    return atmos_bp

def update_atmos(atmocmp, X=1.0, t0=(3.9/100.0), t1=(0.02/100.0), t2=(-0.03/100.0), alpha=-1.7,
                 mol=0.96, BP=782, O3=0.9, H2O=0.9):
    atmocmp.setCoefficients(t0=t0, t1=t1, t2=t2, alpha=alpha, mol=mol, BP=BP, O3=O3, H2O=H2O)
    atmocmp.buildAtmos(secz=X, doplot=True)
    return atmocmp

def atmo_BP(atmocmp):
    atmos_bp = Bandpass(wavelen=atmocmp.wavelen, sb=atmocmp.trans_total)    
    return atmos_bp

###
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

###
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
    print "#Kept %d stars from %s" %(len(starlist), stardir)
    return stars, starlist, temperature, met



if __name__ == "__main__":

    # shift the filters by nothing (standard)
    sys_std = read_hardware(shift_perc=None)
    # shift the filters by one percent
    shift_perc = 0.05
    sys_shift = read_hardware(shift_perc=shift_perc)

    # generate some atmospheres
    atmos = get_atmosDict()

    print atmos.keys()

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
    for atm in atmos.keys():
        mags_std[atm] = {}
        mags_shift[atm] = {}
        for f in filterlist:
            mags_std[atm][f] = numpy.zeros(len(starlist), dtype='float')
            #mags_shift[atm][f] = numpy.zeros(len(starlist), dtype='float')
            i = 0
            for s in starlist:
                mags_std[atm][f][i] = stars[s].calcMag(total_std[atm][f])
                #mags_shift[atm][f][i] = stars[s].calcMag(total_shift[atm][f])
                i = i + 1
    # calculate some colors, in the std bandpass
    gi = mags_std['Standard']['g'] - mags_std['Standard']['i']
    ur = mags_std['Standard']['u'] - mags_std['Standard']['r']
    print gi.min(), gi.max(), ur.min(), ur.max()

    # calculate changes in mag due to
    # changes in the atmosphere
    shifts = {}
    for f in filterlist:
        #shifts[atm][f] = mags_shift[atm][f] - mags_std[atm][f]        
        shifts[f] = mags_std['1.2'][f] - mags_std['2.5'][f]
        # translate to millimags
        shifts[f] = shifts[f] * 1000.0
    print len(gi), len(shifts['u'])


    # make figure
    pylab.figure()
    pylab.subplots_adjust(top=0.93, wspace=0.3, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
    # set colors of data points based on their metallicity (tried first with ur color, but better w/ met)
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
        ax = pylab.subplot(3,2,i)
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
        pylab.ylabel(r"$\Delta$ %s (mmag)" %(f))
        def axis_formatter(x, pos):
            return "%.1f" %(x)
        formatter = pylab.FuncFormatter(axis_formatter)
        ax.yaxis.set_major_formatter(formatter)
        if f == 'u':
            pylab.ylim(-1, 1)
        #elif f=='g':
        #    pylab.ylim(-5, 3)
        #elif f=='i':
        #    pylab.ylim(-1, 1)
        elif f=='y':
            pylab.ylim(-2, 1)
        else:
            pylab.ylim(-1, 1)
        pylab.grid(True)
        i = i + 1
    #pylab.figtext(0.2, 0.95, "Change in magnitude for %s atmosphere and filter shift of %s" %(atmo_choice, shift_perc, "%"))
    pylab.figtext(0.1, 0.95, r"Change in observed magnitude: X=2.0, 10% change in O3/$\tau_0$/$\alpha$, 30% change in H2O")
    #pylab.savefig("delta_mags2.eps", format='eps')

                  
    if False:
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

