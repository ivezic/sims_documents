import os
import numpy
import pylab
import AtmoComp as ac
#import plot_dmags as py
from plot_dmags import *

filters= ('u', 'g', 'r', 'i', 'z', 'y')
colors = ['b', 'g', 'y', 'r', 'm', 'k']

if __name__ == "__main__":

    # Read in default/standard hardware
    sys_std = read_hardware(shift_perc=None)
    # Do you want shifted hardware too? (this is cheap, might as well leave code here)
    shift_perc = -1.0
    sys_shift = read_hardware(shift_perc=shift_perc)

    std_atmos = read_stdatmo()
    # some other atmosphere
    atmocomp, atmos1 = build_atmos(X=1.2, t0=(3.9/100.0), t1=(0.02/100.0), t2=(-0.03/100.0), alpha=-1.7,
                                   mol=0.96, BP=782, O3=1, H2O=1, doPlot=False)
    atmocomp, atmos2= build_atmos(atmocmp=atmocomp, X=1.2, t0=(3.9/100.0), t1=(0.02/100.0),
                                   t2=(-0.03/100.0), alpha=-1.7,
                                   mol=0.96, BP=782, O3=1, H2O=1.5, doPlot=False)
    atmocomp, atmos3 = build_atmos(atmocmp=atmocomp, X=1.2, t0=(3.9/100.0), t1=(0.02/100.0), \
                                   t2=(-0.03/100.0), alpha=-1.7, \
                                   mol=0.96, BP=782, O3=1, H2O=0.5, doPlot=False)

    # plot atmospheres
    pylab.figure()
    pylab.plot(atmos1.wavelen, atmos1.sb, 'k-', label='Nominal H2O')
    pylab.plot(atmos2.wavelen, atmos2.sb, 'r-', label='1.5x H2O')
    pylab.plot(atmos3.wavelen, atmos3.sb, 'g-', label='0.5x H2O')
    pylab.xlim(300, 1100)
    pylab.title('Variable water absorption')
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Transmission')
    leg = pylab.legend(loc='lower right', numpoints=1, fancybox=True, shadow=True)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    
    # combine to create total throughputs
    standard = combine_throughputs(std_atmos, sys_std)

    ##### Water variation in atmosphere
    
    comparison1 = combine_throughputs(atmos1, sys_std)
    comparison2 = combine_throughputs(atmos2, sys_std)
    comparison3 = combine_throughputs(atmos3, sys_std)

    plot_throughputs(comparison1, comparison2, newfig=True, othercolor='r')
    plot_throughputs(comparison1, comparison3, newfig=False, label=False, othercolor='g')
    pylab.xlim(750, 1100)
    pylab.title('Variable water absorption')
    
    # read the kurucz stars
    stars, starlist, temperatures, metallicity, logg = read_kurucz()
    print "Kurucz ", temperatures.max(), temperatures.min(), metallicity.max(), metallicity.min()

    # calculate the standard BP magnitudes for each of these stars
    mags_std = calc_mags(stars, starlist, standard, filterlist)
    gi = calc_stdcolors(mags_std)

    mags1 = calc_adu(stars, starlist, comparison1, filterlist)
    mags2 = calc_adu(stars, starlist, comparison2, filterlist)
    mags3 = calc_adu(stars, starlist, comparison3, filterlist)
    dmags1 = calc_deltamags(mags1, mags2)
    dmags2 = calc_deltamags(mags1, mags3)
    dmags_t = calc_deltamags(mags2, mags3)

    """
    ylims = {}
    xlims = {}
    for f in filterlist:
        ylims[f] = [-10, 0]
        xlims[f] = [-0.06, 0.01]
    """
    ylims = None
    xlims = None

    plot_dmags(gi, dmags_t, [metallicity, logg], 'kurucz', plotfilterlist=('y4',), ylims=ylims,
               titletext="Water variability", newfig=True)

    # do the same for the MLT stars
    do_mlt = True
    if do_mlt:
        mlts, mltlist, mlist, llist, tlist = read_mlt()
        mags_std = calc_mags(mlts, mltlist, standard, filterlist)
        gi = calc_stdcolors(mags_std)
        mags1 = calc_adu(mlts, mltlist, comparison1, filterlist)
        mags2 = calc_adu(mlts, mltlist, comparison2, filterlist)
        mags3 = calc_adu(mlts, mltlist, comparison3, filterlist)
        dmags1 = calc_deltamags(mags1, mags2)
        dmags2 = calc_deltamags(mags1, mags3)
        dmags_t = calc_deltamags(mags2, mags3)
        print "MLTs ", len(mlist), len(llist), len(tlist)
        plot_dmags(gi, dmags_t, [mltlist, mlist, llist, tlist], 'mlt', plotfilterlist=('y4',), ylims=ylims, newfig=False)

    do_sn = True
    if do_sn:
        sns, snlist, days, redshifts = read_sn()
        mags_std = calc_mags(sns, snlist, standard, filterlist)
        gi = calc_stdcolors(mags_std)
        mags1 = calc_adu(sns, snlist, comparison1, filterlist)
        mags2 = calc_adu(sns, snlist, comparison2, filterlist)
        mags3 = calc_adu(sns, snlist, comparison3, filterlist)
        dmags1 = calc_deltamags(mags1, mags2)
        dmags2 = calc_deltamags(mags1, mags3)
        dmags_t = calc_deltamags(mags2, mags3)
        print "SN " , snlist
        plot_dmags(gi, dmags_t, [snlist,], 'sn', plotfilterlist=('y4',), ylims=ylims, newfig=False)


    ####  BANDPASS SHIFT
    comparison1 = combine_throughputs(atmos1, sys_std)
    comparison2 = combine_throughputs(atmos1, sys_shift)
    filters1 = read_filtersonly(shift_perc=None)
    filters2 = read_filtersonly(shift_perc=shift_perc)
    # plot filter bandpasses
    pylab.figure()
    pylab.plot(filters1['z'].wavelen, filters1['z'].sb, 'k-', label='Nominal BP')
    pylab.plot(filters2['z'].wavelen, filters2['z'].sb, 'r-', label='1% Shift')
    pylab.plot(filters1['y4'].wavelen, filters1['y4'].sb, 'k-')
    pylab.plot(filters2['y4'].wavelen, filters2['y4'].sb, 'r-')
    pylab.xlim(750, 1100)
    pylab.title('Filter Bandpass Shift')
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Transmission')
    leg = pylab.legend(loc='lower right', numpoints=1, fancybox=True, shadow=True)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
        

    plot_throughputs(comparison1, comparison2, newfig=True, othercolor='r')
    pylab.xlim(750, 1100)
    pylab.title('Filter Bandpass Shift')

    # read the kurucz stars
    #stars, starlist, temperatures, metallicity, logg = read_kurucz()
    print "Kurucz ", temperatures.max(), temperatures.min(), metallicity.max(), metallicity.min()

    # calculate the standard BP magnitudes for each of these stars
    mags_std = calc_mags(stars, starlist, standard, filterlist)
    gi = calc_stdcolors(mags_std)

    mags1 = calc_adu(stars, starlist, comparison1, filterlist)
    mags2 = calc_adu(stars, starlist, comparison2, filterlist)
    dmags_t = calc_deltamags(mags1, mags2)

    """
    ylims = {}
    xlims = {}
    for f in filterlist:
        ylims[f] = [-10, 0]
        xlims[f] = [-0.06, 0.01]
    """
    ylims = None
    xlims = None

    plot_dmags(gi, dmags_t, [metallicity, logg], 'kurucz', plotfilterlist=('y4',), ylims=ylims,
               titletext="Bandpass Shift", newfig=True)

    # do the same for the MLT stars
    do_mlt = True
    if do_mlt:
        # mlts, mltlist, mlist, llist, tlist = read_mlt()
        mags_std = calc_mags(mlts, mltlist, standard, filterlist)
        gi = calc_stdcolors(mags_std)
        mags1 = calc_adu(mlts, mltlist, comparison1, filterlist)
        mags2 = calc_adu(mlts, mltlist, comparison2, filterlist)
        dmags_t = calc_deltamags(mags1, mags2)
        print "MLTs ", len(mlist), len(llist), len(tlist)
        plot_dmags(gi, dmags_t, [mltlist, mlist, llist, tlist], 'mlt', plotfilterlist=('y4',), ylims=ylims, newfig=False)

    do_sn = True
    if do_sn:
        #sns, snlist, days, redshifts = read_sn()
        mags_std = calc_mags(sns, snlist, standard, filterlist)
        gi = calc_stdcolors(mags_std)
        mags1 = calc_adu(sns, snlist, comparison1, filterlist)
        mags2 = calc_adu(sns, snlist, comparison2, filterlist)
        dmags_t = calc_deltamags(mags1, mags2)
        print "SN " , snlist
        plot_dmags(gi, dmags_t, [snlist,], 'sn', plotfilterlist=('y4',), ylims=ylims, newfig=False)



    pylab.show()
    
