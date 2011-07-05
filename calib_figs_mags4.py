import os
import numpy
import pylab
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import lsst.sims.catalogs.measures.photometry.photUtils as photUtils
import AtmoComp as ac

WMIN = 300
WMAX = 1200
WSTEP = 0.1
filterlist = ('u', 'g', 'r', 'i', 'z', 'y3')
colors = ['b', 'g', 'y', 'r', 'm', 'k']
def color_counter_next(i):
    i = i + 1
    if i == len(colors):
        i = 0
    return i

def read_hardware(shift_perc=None):
    # read system (hardware) transmission, return dictionary of system hardware (keyed to filter)
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
    # read the filter throughput curves only (called from read_hardware as well)
    # apply a shift of +shift_perc/100 * eff_wavelength to the wavelengths of the filter.
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
    # read (or build) a series of atmosphere files, keyed with desired names for the atmosphere.
    atmos = {}
    # First read the normal, standardized atmosphere from a file.
    atmos['Standard'] = read_stdatmo()
    # Then set up the atmospheres that should be build from templates. Initialize/read the templates.
    atmocmp = ac.AtmoComp()    
    for key in ('1', '2'):
        X = 1.2  # float(key)
        if key=='1':
            atmocmp, atmos[key] = get_atmos(atmocmp, X=X, t0=5.6/100.0, alpha=-1.8, O3=1.5, H2O=1.3)
        if key=='2':
            atmocmp, atmos[key] = get_atmos(atmocmp, X=X, t0=0.2/100.0, alpha=-0.5, O3=0.6, H2O=0.5)
        # examples
        # max atmo = get_atmos(atmocmp, X=X, t0=5.6/100.0, alpha=-1.8, O3=1.5, H2O=1.3)
        # min atmo = get_atmos(atmocmp, X=X, t0=0.2/100.0, alpha=-0.5, O3=0.6, H2O=0.5)
        # 30p atmo = get_atmos(atmocmp, X=2.5, t0=(0.8/100), alpha=-1.0, O3=0.9, H2O=0.8)
        # 30p atmo = get_atmos(atmocmp, X=2.5, t0=(2.4/100.0), alpha=-1.4, O3=1.17, H2O=1.04)
        # 10p/30 atmo = get_atmos(atmocmp, X=X, t0=(0.8/100.0), alpha=-1.0, O3=0.9, H2O=0.8)
        # 10p/30 atmo = get_atmos(atmocmp, X=X, t0=(1.3/100.0), alpha=-1.13, O3=0.99, H2O=1.04)
    return atmos

def read_stdatmo():
    atmosdir = "."
    atmos_bp = Bandpass()
    atmos_bp.readThroughput(os.path.join(atmosdir, "atmos_std.dat"))
    return atmos_bp

def get_atmos(atmocmp, X=1.0, t0=(3.9/100.0), t1=(0.02/100.0), t2=(-0.03/100.0), alpha=-1.7,
              mol=0.96, BP=782, O3=0.9, H2O=0.9):
    atmocmp.setCoefficients(t0=t0, t1=t1, t2=t2, alpha=alpha, mol=mol, BP=BP, O3=O3, H2O=H2O)
    atmocmp.buildAtmos(secz=X, doplot=True)
    atmos_bp = Bandpass(wavelen=atmocmp.wavelen, sb=atmocmp.trans_total)    
    return atmocmp, atmos_bp

###
def combine_throughputs(atmos, sys):
    # Set up the total throughput for this system bandpass, using the variety of atmospheres.
    # returns a dictionary of total throughputs, keyed for each key in atmosphere and then by filter.
    total = {}
    for key in atmos.keys():
        total[key] = {}
        for f in filterlist:
            wavelen, sb = sys[f].multiplyThroughputs(atmos[key].wavelen, atmos[key].sb)
            total[key][f] = Bandpass(wavelen, sb)
            total[key][f].sbTophi()
    return total


def plot_throughputs(bpDict1, bpDict2):
    pylab.figure()
    i = 0
    for f in filterlist:
        pylab.plot(bpDict1[f].wavelen, bpDict1[f].sb, colors[i]+"-", label=f)
        pylab.plot(bpDict2[f].wavelen, bpDict2[f].sb, colors[i]+":")
        i = color_counter_next(i)
    pylab.ylim(0, 0.8)
    pylab.xlim(300, 1200)
    pylab.xlabel("Wavelength (nm)")
    pylab.ylabel("Transmission")
    pylab.grid(True)
    pylab.legend(numpoints=1, fancybox=True, loc='upper right')
    return

###
def read_kurucz():
    # read kurucz model MS, g40 stars SEDs
    homedir = os.getenv("HOME")    
    stardir = os.path.join(homedir, "seds/kurucz_r")
    allfilelist = os.listdir(stardir)
    starlist = []
    # make preliminary cut for ms, g40 stars
    for filename in allfilelist:
        if filename[-3:] == 'g40':
            starlist.append(filename)
    atemperature = []
    amet = []
    starlist2 = []
    # make secondary cut for stars with temperature > 5000 deg
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
    # actually read the stars SEDS from disk
    stars = {}
    for s in starlist:
        stars[s] = Sed()
        stars[s].readSED_flambda(os.path.join(stardir,s))
    print "# Read %d MS stars from %s" %(len(starlist), stardir)
    # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
    for s in starlist:
        stars[s].synchronizeSED(wavelen_min=WMIN, wavelen_max=WMAX, wavelen_step=WSTEP)
    return stars, starlist, temperature, met

def read_mlt():
    # read mlt stars - only keep 'm's
    # find the filenames and mark 'm', 'l', 't' stars separately
    homedir = os.getenv("HOME")
    mltdir = os.path.join(homedir, "seds/mlt")
    allfilelist = os.listdir(mltdir)
    mltlist = []
    mlist = []
    llist = []
    tlist = []
    for filename in allfilelist:
        if filename.endswith('.dat') & filename.startswith('m'):
            mlist.append(filename)
        elif filename.endswith('.dat') & filename.startswith('L'):
            llist.append(filename)
        elif filename.startswith('burrows'):
            tlist.append(filename)
    mltlist = mlist # + llist + tlist
    # read the mlt seds from disk
    mlts = {}
    for s in mltlist:
        mlts[s] = Sed()
        mlts[s].readSED_flambda(os.path.join(mltdir, s))
    print "# Read %d mlt stars from %s" %(len(mltlist), mltdir)
    # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
    for s in mltlist:
        mlts[s].synchronizeSED(wavelen_min=WMIN, wavelen_max=WMAX, wavelen_step=WSTEP)
    return mlts, mltlist, mlist, llist, tlist

def read_quasar():
    # read quasar spectra and redshift
    homedir = os.getenv("HOME")
    quasardir = os.path.join(homedir, "seds/quasar")
    # read zero redshift quasar
    base = Sed()
    base.readSED_flambda(os.path.join(quasardir, "quasar.dat"))
    # redshift 
    #redshifts = [0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5]
    #redshifts = numpy.array(redshifts)
    redshifts= numpy.arange(0, 2.8, 0.1)
    quasars = {}
    for z in redshifts:
        wavelen, flambda = base.redshiftSED(z, wavelen=base.wavelen, flambda=base.flambda)
        quasars[z] = Sed(wavelen=wavelen, flambda=flambda)
    print "# Generated %d quasars at redshifts between %f and %f" %(len(redshifts), redshifts.min(), redshifts.max())
    # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
    for z in redshifts:
        quasars[z].synchronizeSED(wavelen_min=WMIN, wavelen_max=WMAX, wavelen_step=WSTEP)
    return quasars, redshifts

def read_whitedwarf():
    # read white dwarf bergeron models
    homedir = os.getenv("HOME")
    whitedwarfdir = os.path.join(homedir, "seds/white_dwarfs_r")
    # get the H dwarfs
    Hdir = os.path.join(whitedwarfdir, "H")
    allfilelist = os.listdir(Hdir)
    hlist = []
    temperatures = []
    loggs = []
    for filename in allfilelist:
        if filename.startswith('bergeron'):
            tmp = filename.split('_')
            temperature = float(tmp[1])
            logg = float(tmp[2].split('.')[0])
            logg = logg/10.0
            if (logg > 7.0) & (temperature>5000):
                hlist.append(filename)
                temperatures.append(temperature)
                loggs.append(logg)
    Hedir = os.path.join(whitedwarfdir, "He")
    allfilelist = os.listdir(Hedir)
    helist = []
    for filename in allfilelist:
        if filename.startswith('bergeron_He'):
            tmp = filename.split('_')
            temperature = float(tmp[2])
            logg = float(tmp[3].split('.')[0])
            logg = logg/10.0
            if (logg > 7.0) & (temperature>5000):
                helist.append(filename)                        
                temperatures.append(temperature)
                loggs.append(logg)
    temperatures = numpy.array(temperatures)
    loggs = numpy.array(loggs)
    wdlist = hlist + helist
    wds = {}
    for w in wdlist:
        wds[w] = Sed()
        if w in hlist:
            wds[w].readSED_flambda(os.path.join(Hdir, w))
        if w in helist:
            wds[w].readSED_flambda(os.path.join(Hedir, w))
    # synchronize seds for faster mag calcs later
    for w in wdlist:
        wds[w].synchronizeSED(wavelen_min=WMIN, wavelen_max=WMAX, wavelen_step=WSTEP)
    return wds, wdlist, hlist, helist, temperatures, loggs

def read_sn():
    # read sn spectra and redshift
    homedir = os.getenv("HOME")
    sndir = os.path.join(homedir, "seds/sn")
    allfilelist = os.listdir(sndir)
    snlist = []
    days = ['0', '20', '40']
    #redshifts = [0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5]
    #redshifts = numpy.array(redshifts)
    redshifts= numpy.arange(0, 1.0, 0.1)
    # pull out the filenames we want
    for filename in allfilelist:
        if filename.endswith('.dat') & filename.startswith('sn1a_'):
            #snlist.append(filename)
            tmp = filename.split('_')
            day = tmp[1].split('.')[0]
            if day in days:
                snlist.append(filename)
    # read base SEDs for these days
    sns_base = {}
    for r in zip(snlist, days):
        day = r[1]
        sns_base[day] = Sed()
        sns_base[day].readSED_flambda(os.path.join(sndir, r[0]))
    # and redshift
    sns = {}
    snlist = []
    for d in days:        
        for z in redshifts:
            sn_name = "%d_%.1f" %(int(d), z)
            wavelen, flambda = sns_base[d].redshiftSED(z, wavelen=sns_base[d].wavelen, flambda=sns_base[d].flambda)
            sns[sn_name] = Sed(wavelen=wavelen, flambda=flambda)
            snlist.append(sn_name)
    print "# Generated %d sn's at redshifts between %f and %f on days %s" %(len(snlist),
                                                                            redshifts.min(), redshifts.max(), days)
    # resample onto the standard bandpass for Bandpass obj's and calculate fnu to speed later calculations
    for s in snlist:
        sns[s].synchronizeSED(wavelen_min=WMIN, wavelen_max=WMAX, wavelen_step=WSTEP)
    return sns, snlist, days, redshifts

###

def calc_mags(seds, sedkeylist, bpDict):
    # calculate magnitudes for all sed objects using bpDict (e.g. total['Standard'] or shifted['1'])
    # pass the sedkeylist so you know what order the magnitudes are arranged in
    mags = {}
    for f in filterlist:
        mags[f] = numpy.zeros(len(sedkeylist), dtype='float')
        i = 0
        for key in sedkeylist:
            mags[f][i] = seds[key].calcMag(bpDict[f])
            i = i + 1
    return mags

def calc_stdcolors(mags_std):
    # calculate some colors in the standard atmosphere, should be also standard bandpass, not shifted)
    gi = mags_std['g'] - mags_std['i']
    return gi

def plot_colorcolor(seds, sedkeylist, sedcolorkey, sedtype, bpDict, titletext=None, newfig=True):
    if newfig:
        pylab.figure()
    mags = calc_mags(seds, sedkeylist, bpDict)
    gi = calc_stdcolors(mags)
    magcolors = {}
    colorlabels = []
    colorslabels[0] = None
    for i in range(1-7):
        colorlabels[i] = filterlist[i-1] + '-' + filterlist[i]
        magscolors[colorlabels[i]] = mags[filterlist[i-1]] - mags[filterlist[i]]
    # colors = ug, gr, ri, iz, zy
    if sedtype == 'kurucz':
        print "# plotting kurucz"
       # set colors of data points based on their metallicity
        metallicity = sedcolorkey
        metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
        metbinsize = abs(sedcolorkey.min() - sedcolorkey.max())/6.0
        metbins = numpy.arange(sedcolorkey.min(), sedcolorkey.max() + metbinsize, metbinsize)
        print metbinsize, metbins
        i = 1
        # use a different subplot for each color/color combo
        for i in range(len(colorlabels-1)):
            ax = pylab.subplot(3,2,i)
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
                mcolor = metcolors[metidx]
                pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], mcolor+'.')
            i = i + 1
        ax = pylab.subplot(3,2,7)
        for metidx in range(len(metbins)):
            condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
            mcolor = metcolors[metidx]
            pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], mcolor+'.')
    elif sedtype == 'quasar':
        print "# looking at quasars"
        # set the colors of data points based on their redshift
        redshift = sedcolorkey
        redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
        redbinsize = 0.5
        redbins = numpy.arange(0.0, 3.0+redbinsize, redbinsize)
        print redbinsize, redbins, redcolors
        i = 1
        # for each filter, use a different subplot
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for redidx in range(len(redbins)):
                condition =((redshift>=redbins[redidx]) & (redshift<=redbins[redidx]+redbinsize))
                rcolor = redcolors[redidx]
                pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], rcolor+'o')
            i = i + 1
        ax = pylab.subplot(3,2,7)
        for redidx in range(len(redbins)):
            condition =((redshift>=redbins[redidx]) & (redshift<=redbins[redidx]+redbinsize))
            rcolor = redcolors[redidx]
            pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], rcolor+'o')
    elif sedtype == 'mlt':
        print "# looking at MLT"
        # set the colors of data points based on their type
        mltlist = sedcolorkey[0]
        mlist = sedcolorkey[1]
        llist = sedcolorkey[2]
        tlist = sedcolorkey[3]
        i = 1
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for j in range(len(mltlist)):
                if (mltlist[j] in mlist):
                    pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], 'bx')
                elif (mltlist[j] in llist):
                    pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], 'gx')
                elif (mltlist[j] in tlist):
                    pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], 'mx')
            i = i + 1
        ax = pylab.subplot(3,2,7)
        for j in range(len(mltlist)):
            if (mltlist[j] in mlist):
                pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], 'bx')
            elif (mltlist[j] in llist):
                pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], 'gx')
            elif (mltlist[j] in tlist):
                pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], 'mx')
    elif sedtype == 'white_dwarf':
        print "# looking at White Dwarf"
        # set the colors of data points based on their type
        wdlist = sedcolorkey[0]
        hlist = sedcolorkey[1]
        helist = sedcolorkey[2]
        i = 1
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for j in range(len(wdlist)):
                if (wdlist[j] in hlist):
                    pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], 'y+')
                elif (wdlist[j] in helist):
                    pylab.plot(magscolors[colorlabels[i]][condition], magscolors[colorlabels[i+1]][f][condition], 'y+')
            i = i + 1
        ax = pylab.subplot(3,2,7)
        for j in range(len(wdlist)):
            if (wdlist[j] in hlist):
                pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], 'y+')
            elif (wdlist[j] in helist):
                pylab.plot(gi[condition], magscolors[colorlabels[i+1]][f][condition], 'y+')
    # set up generic items
    for i in range(1, 7):
        f = filterlist[i-1]
        ax = pylab.subplot(3,2,i)
        #pylab.xlabel("g-i")
        #pylab.ylabel(r"$\Delta$ %s (mmag)" %(f))
        def axis_formatter(x, pos):
            return "%.1f" %(x)
        formatter = pylab.FuncFormatter(axis_formatter)
        ax.yaxis.set_major_formatter(formatter)
        # set axes limits
        if ylims == None:
            pass
        else:
            try:
                pylab.ylim(ylims[f][0], ylims[f][1])
            except KeyError:
                pass
    # put a grid in the background
    if newfig:
        for i in range(1, 7):
            ax = pylab.subplot(3, 2, i)
            pylab.grid(True)
    #pylab.suptitle(titletext)
    #pylab.savefig("delta_mags2.eps", format='eps')
    return

                
def calc_deltamags(seds, sedkeylist, bpDict1, bpDict2):
    # calculate change in magnitude due to change in some parameter difference between bpDict1/2
    # give this function the bp dictionary stripped down to the filterlist level
    #  (e.g. bpDict1 =total['Standard'], bpDict2 = total['1'])
    dmags = {}
    for f in filterlist:
        dmags[f] = numpy.zeros(len(sedkeylist), dtype='float')
        i = 0
        for key in sedkeylist:
            # calculate difference in magnitude
            dmags[f][i] = (seds[key].calcMag(bpDict1[f]) - seds[key].calcMag(bpDict2[f])) * 1000.0
            # converted to mmag
            i = i + 1
    return dmags

def plot_dmags(gi, dmags, sedcolorkey, sedtype, titletext=None, ylims=None, xlims=None, newfig=True):
    # ylims = dictionary (ylim['u'] = [0, 0])
    sedtypes = ['kurucz', 'mlt', 'quasar', 'white_dwarf', 'sn']
    # make figure of change in magnitude
    if newfig:
        pylab.figure()
    pylab.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
    # For Kurucz models
    if sedtype == 'kurucz':
        print "# looking at kurucz"
        # set colors of data points based on their metallicity
        metallicity = sedcolorkey
        metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
        metbinsize = abs(sedcolorkey.min() - sedcolorkey.max())/6.0
        metbins = numpy.arange(sedcolorkey.min(), sedcolorkey.max() + metbinsize, metbinsize)
        print metbinsize, metbins
        i = 1
        # for each filter, use a different subplot
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
                mcolor = metcolors[metidx]
                pylab.plot(gi[condition], dmags[f][condition], mcolor+'.')
            i = i + 1
    elif sedtype == 'quasar':
        print "# looking at quasars"
        # set the colors of data points based on their redshift
        redshift = sedcolorkey
        redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
        redbinsize = 0.5
        redbins = numpy.arange(0.0, 3.0+redbinsize, redbinsize)
        print redbinsize, redbins, redcolors
        i = 1
        # for each filter, use a different subplot
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for redidx in range(len(redbins)):
                condition =((redshift>=redbins[redidx]) & (redshift<=redbins[redidx]+redbinsize))
                rcolor = redcolors[redidx]
                pylab.plot(gi[condition], dmags[f][condition], rcolor+'o')
            i = i + 1
    elif sedtype == 'mlt':
        print "# looking at MLT"
        # set the colors of data points based on their type
        mltlist = sedcolorkey[0]
        mlist = sedcolorkey[1]
        llist = sedcolorkey[2]
        tlist = sedcolorkey[3]
        i = 1
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for j in range(len(mltlist)):
                if (mltlist[j] in mlist):
                    pylab.plot(gi[j], dmags[f][j], 'bx')
                elif (mltlist[j] in llist):
                    pylab.plot(gi[j], dmags[f][j], 'gx')
                elif (mltlist[j] in tlist):
                    pylab.plot(gi[j], dmags[f][j], 'mx')
            i = i + 1
    elif sedtype == 'white_dwarf':
        print "# looking at White Dwarf"
        # set the colors of data points based on their type
        wdlist = sedcolorkey[0]
        hlist = sedcolorkey[1]
        helist = sedcolorkey[2]
        i = 1
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            for j in range(len(wdlist)):
                if (wdlist[j] in hlist):
                    pylab.plot(gi[j], dmags[f][j], 'y+')
                elif (wdlist[j] in helist):
                    pylab.plot(gi[j], dmags[f][j], 'y+')
            i = i + 1
    elif sedtype == 'sn':
        print "# looking at SN"
        # set the color of data points based on their redshift and day
        snlist = sedcolorkey
        redcolors = ['b', 'b', 'g', 'g', 'r', 'r' ,'m', 'm']
        redbinsize = 0.18
        redbins = numpy.arange(0.0, 1.0+redbinsize, redbinsize)
        print redbinsize, redbins, redcolors
        day_symbol = {'0':'s', '20':'s', '40':'s'}
        i = 1
        for f in filterlist:
            ax = pylab.subplot(3,2,i)
            j = 0
            for s in snlist:
                day, redshift = s.split('_')
                redshift = float(redshift)
                redidx = int(redshift / redbinsize)
                pylab.plot(gi[j], dmags[f][j], redcolors[redidx]+day_symbol[day])
                j = j + 1
            i = i + 1
    # set up generic items
    for i in range(1, 7):
        f = filterlist[i-1]
        ax = pylab.subplot(3,2,i)
        pylab.xlabel("g-i")
        pylab.ylabel(r"$\Delta$ %s (mmag)" %(f))
        def axis_formatter(x, pos):
            return "%.1f" %(x)
        formatter = pylab.FuncFormatter(axis_formatter)
        ax.yaxis.set_major_formatter(formatter)
        # set axes limits
        if ylims == None:
            pass
        else:
            try:
                pylab.ylim(ylims[f][0], ylims[f][1])
            except KeyError:
                pass
        if xlims == None:
            pass
        else:
            try:
                pylab.xlim(xlims[f][0], xlims[f][1])
            except KeyError:
                pass
    # put a grid in the background
    if newfig:
        for i in range(1, 7):
            ax = pylab.subplot(3, 2, i)
            pylab.grid(True)            
        pylab.suptitle(titletext)
    #pylab.savefig("delta_mags2.eps", format='eps')
    return

if __name__ == "__main__":

    # Read in standard hardware
    sys_std = read_hardware(shift_perc=None)
    # shift the filters by some percent
    shift_perc = 1.0
    sys_shift = read_hardware(shift_perc=shift_perc)

    # generate some atmospheres
    atmos = get_atmosDict()
    print atmos.keys()

    # combine to create total throughputs 
    total_std = combine_throughputs(atmos, sys_std)    
    total_shift = combine_throughputs(atmos, sys_shift)

    # read the kurucz stars
    stars, starlist, temperatures, metallicity = read_kurucz()
    print "Kurucz ", temperatures.max(), temperatures.min(), metallicity.max(), metallicity.min()

    # calculate the standard BP magnitudes for each of these stars
    mags_std = calc_mags(stars, starlist, total_std['Standard'])
    gi = calc_stdcolors(mags_std)

    # now calculate change in magnitude between the two quantities we want to compare
    t1 = total_std['1']
    t2 = total_shift['1']
    #titletext = "Maximum changes in atmosphere, all SED types"
    #titletext = "10% / 30% H2O atmosphere changes"
    titletext = "1% filter shift"
    dmags = calc_deltamags(stars, starlist, t1, t2)
    for f in filterlist:        
        print f, dmags[f].min(), dmags[f].max()

    #ylims['u'] = [-100, 20]
    #ylims['g'] = [-40, 20]
    #ylims['r'] = [-20, 25]
    #ylims['i'] = [-20, 25]
    #ylims['z'] = [-20, 25]
    #ylims['y'] = [-20, 25]
    ylims = {}
    xlims = {}
    for f in filterlist:
        ylims[f] = [-10, 0]
        xlims[f] = [-0.06, 0.01]

    #plot_throughputs(t1, t2)

    plot_dmags(gi, dmags, metallicity, 'kurucz', titletext=titletext, ylims=ylims, newfig=True)

    #pylab.show()
    #exit()

    # do the same for the MLT stars
    do_mlt = True
    if do_mlt:
        mlts, mltlist, mlist, llist, tlist = read_mlt()
        mags_std = calc_mags(mlts, mltlist, total_std['Standard'])
        gi = calc_stdcolors(mags_std)
        dmags = calc_deltamags(mlts, mltlist, t1, t2)
        print "MLTs ", len(mlist), len(llist), len(tlist)
        for f in filterlist:        
            print f, dmags[f].min(), dmags[f].max()
        plot_dmags(gi, dmags, [mltlist, mlist, llist, tlist], 'mlt', titletext=titletext, ylims=ylims, newfig=False)

    # and for the quasars
    do_quasar = True
    if do_quasar:
        quasars, redshifts = read_quasar()
        mags_std = calc_mags(quasars, redshifts, total_std['Standard'])
        gi  = calc_stdcolors(mags_std)
        dmags = calc_deltamags(quasars, redshifts, t1, t2)
        print "Quasars ", redshifts
        for f in filterlist:
            print f, dmags[f].min(), dmags[f].max()
        plot_dmags(gi, dmags, redshifts, 'quasar', titletext=titletext, ylims=ylims, newfig=False)

    do_whitedwarf = False
    if do_whitedwarf:
        wds, wdlist, hlist, helist, temperatures, loggs = read_whitedwarf()
        mags_std = calc_mags(wds, wdlist, total_std['Standard'])
        gi  = calc_stdcolors(mags_std)
        dmags = calc_deltamags(wds, wdlist, t1, t2)
        print "White dwarfs ", len(hlist), len(helist)
        for f in filterlist:
            print f, dmags[f].min(), dmags[f].max()
        plot_dmags(gi, dmags, [wdlist, hlist, helist], 'white_dwarf', titletext=titletext, ylims=ylims, newfig=False)


    do_sn = True
    if do_sn:
        sns, snlist, days, redshifts = read_sn()
        mags_std = calc_mags(sns, snlist, total_std['Standard'])
        gi = calc_stdcolors(mags_std)
        dmags = calc_deltamags(sns, snlist, t1, t2)
        print "SN " , snlist
        for f in filterlist:
            print f, dmags[f].min(), dmags[f].max()
        plot_dmags(gi, dmags, snlist, 'sn', titletext=titletext, ylims=ylims, xlims=xlims, newfig=False)

        
    pylab.show()

