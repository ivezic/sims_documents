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

def read_hardware():
    # read system (hardware) transmission)
    filterdir = os.getenv("LSST_THROUGHPUTS_DEFAULT")
    hardware = ("detector.dat", "m1.dat", "m2.dat", "m3.dat", "lens1.dat", "lens2.dat", "lens3.dat")
    sys_std = {}
    # Read in all the standard components
    for f in filterlist:
        sys_std[f] = Bandpass()
        tlist = []
        for t in hardware:
            tlist.append(os.path.join(filterdir, t))
        tlist.append(os.path.join(filterdir, "filter_" + f + ".dat"))
        sys_std[f].readThroughputList(tlist)
    sys_edge = {}
    # Read in the standard components, except shift the filter by 1% of the central wavelength near the edge
    for f in filterlist:
        sys_edge[f] = Bandpass()
        tlist = []
        for t in hardware:
            tlist.append(os.path.join(filterdir, t))
        sys_edge[f].readThroughputList(tlist)
        tmpfilter = Bandpass()
        tmpfilter.readThroughput(os.path.join(filterdir, "filter_" + f + ".dat"))
        effwavelenphi, effwavelensb = tmpfilter.calcEffWavelen()
        shift = effwavelensb * 0.01
        tmpfilter.wavelen = tmpfilter.wavelen + shift
        tmpfilter.resampleBandpass()
        sys_edge[f].wavelen, sys_edge[f].sb = sys_edge[f].multiplyThroughputs(tmpfilter.wavelen, tmpfilter.sb)
    for f in filterlist:
        junk, effwavelen = sys_std[f].calcEffWavelen()
        junk, effwavelen_edge = sys_edge[f].calcEffWavelen()
        print f, effwavelen, effwavelen_edge, effwavelen/effwavelen_edge
    return sys_std, sys_edge

def read_atmos():
    # read atmosphere throughputs, make plots
    atmosdir = "."
    atmos = {}
    key = "Standard"
    atmos[key]= Bandpass()
    atmos[key].readThroughput(os.path.join(atmosdir, "atmos_std.dat"))
    key = "X=1.0, H2O=0.7"
    atmos[key] = Bandpass()
    atmos[key].readThroughput(os.path.join(atmosdir, "atmos_85233890.dat"))
    key = "X=1.5, H2O=1.4"
    atmos[key] = Bandpass()
    atmos[key].readThroughput(os.path.join(atmosdir, "atmos_85488653.dat"))
    return atmos

def combine_throughputs(atmos, sys_std, sys_edge):
    # Set up all the 'total' throughputs
    total = {}
    key = 'standard'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys_std[f].multiplyThroughputs(atmos['Standard'].wavelen, atmos['Standard'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'standard, edge'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys_edge[f].multiplyThroughputs(atmos['Standard'].wavelen, atmos['Standard'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'low X, center'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys_std[f].multiplyThroughputs(atmos['X=1.0, H2O=0.7'].wavelen, atmos['X=1.0, H2O=0.7'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'hi X, center'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys_std[f].multiplyThroughputs(atmos['X=1.5, H2O=1.4'].wavelen, atmos['X=1.5, H2O=1.4'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'low X, edge'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys_edge[f].multiplyThroughputs(atmos['X=1.0, H2O=0.7'].wavelen, atmos['X=1.0, H2O=0.7'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    key = 'hi X, edge'
    total[key] = {}
    for f in filterlist:
        wavelen, sb = sys_edge[f].multiplyThroughputs(atmos['X=1.5, H2O=1.4'].wavelen, atmos['X=1.5, H2O=1.4'].sb)
        total[key][f] = Bandpass(wavelen, sb)
        total[key][f].sbTophi()
    pylab.figure()
    for key in total.keys():
        i = 0
        for f in filterlist:
            pylab.plot(total[key][f].wavelen, total[key][f].phi, colors[i])
            i = color_counter_next(i)
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


def calc_mags(stars, total):
    tkey = "standard"
    mags_std = {}
    mags_std['blue'] = {}
    mags_std['red'] = {}
    writestring1  = "Standard atm, std sys  &  red "
    writestring2 =  "Standard atm, std sys  &  blue "
    for f in filterlist:
        for skey in stars.keys():
            mags_std[skey][f] = stars[skey].calcMag(total[tkey][f])
        writestring1 = writestring1 + "& %.3f " %(mags_std['red'][f])
        writestring2 = writestring2 + "& %.3f " %(mags_std['blue'][f])
    print writestring1 + "\\\\"
    print writestring2 + "\\\\ \\hline \\hline"
    tkey = "standard, edge"
    mags_delta = {}
    mags_delta['blue'] = {}
    mags_delta['red'] = {}
    key = "std atm, +1% sys"
    mags_delta['blue'][key] = {}
    mags_delta['red'][key] = {}
    writestring1  = "Standard atm, +1\% sys shift & red  "
    writestring2  = "Standard atm, +1\% sys shift & blue  "
    for f in filterlist:
        for skey in stars.keys():
            mags_delta[skey][key][f] = stars[skey].calcMag(total[tkey][f]) - mags_std[skey][f]
        writestring1 = writestring1 + "& %.3f " %(mags_delta['red'][key][f])
        writestring2 = writestring2 + "& %.3f " %(mags_delta['blue'][key][f])
    print writestring1 + "\\\\"
    print writestring2 + "\\\\ \\hline"
    tkey = "low X, center"
    key = "X=1.0, std sys"
    mags_delta['blue'][key] = {}
    mags_delta['red'][key] = {}
    writestring1  = "X=1.0, std sys & red  "
    writestring2  = "X=1.0, std sys & blue  "
    for f in filterlist:
        for skey in stars.keys():
            mags_delta[skey][key][f] = stars[skey].calcMag(total[tkey][f]) - mags_std[skey][f]
        writestring1 = writestring1 + "& %.3f " %(mags_delta['red'][key][f])
        writestring2 = writestring2 + "& %.3f " %(mags_delta['blue'][key][f])
    print writestring1 + "\\\\"
    print writestring2 + "\\\\ \\hline"
    tkey = "low X, edge"
    key = "X=1.0, +1% sys"
    mags_delta['blue'][key] = {}
    mags_delta['red'][key] = {}
    writestring1  = "X=1.0, +1\% sys shift & red "
    writestring2  = "X=1.0, +1\% sys shift & blue "
    for f in filterlist:
        for skey in stars.keys():
            mags_delta[skey][key][f] = stars[skey].calcMag(total[tkey][f]) - mags_std[skey][f]
        writestring1 = writestring1 + "& %.3f " %(mags_delta['red'][key][f])
        writestring2 = writestring2 + "& %.3f " %(mags_delta['blue'][key][f])
    print writestring1 + "\\\\"
    print writestring2 + "\\\\ \\hline"
    tkey = "hi X, center"
    key = "X=1.5, std sys"
    mags_delta['blue'][key] = {}
    mags_delta['red'][key] = {}
    writestring1  = "X=1.5, std sys  & red "
    writestring2  = "X=1.5, std sys  & blue "
    for f in filterlist:
        for skey in stars.keys():
            mags_delta[skey][key][f] = stars[skey].calcMag(total[tkey][f]) - mags_std[skey][f]
        writestring1 = writestring1 + "& %.3f " %(mags_delta['red'][key][f])
        writestring2 = writestring2 + "& %.3f " %(mags_delta['blue'][key][f])
    print writestring1 + "\\\\"
    print writestring2 + "\\\\ \\hline"
    tkey = "hi X, edge"
    key = "X=1.5, +1% sys"
    mags_delta['blue'][key] = {}
    mags_delta['red'][key] = {}
    writestring1 = "X=1.5, +1\% sys shift & red "
    writestring2  = "X=1.5, +1\% sys shift & blue "
    for f in filterlist:
        for skey in stars.keys():
            mags_delta[skey][key][f] = stars[skey].calcMag(total[tkey][f]) - mags_std[skey][f]
        writestring1 = writestring1 + "& %.3f " %(mags_delta['red'][key][f])
        writestring2 = writestring2 + "& %.3f " %(mags_delta['blue'][key][f])
    print writestring1 + "\\\\"
    print writestring2 + "\\\\ \\hline"
    return mags_delta

def make_figures(stars, sys_std, sys_edge, atmos, total, mags_delta):
    pylab.figure()
    pylab.subplots_adjust(top=0.93, wspace=0.15, hspace=0.3, bottom=0.08, left=0.09, right=0.91)
    # make plot of flux for these stars
    ax = pylab.subplot(532)
    scolor = ['b', 'r']
    i =0 
    for key in stars.keys():
        ax.plot(stars[key].wavelen, stars[key].flambda*1e15, scolor[i], label=key)
        i = i +1
    #pylab.legend(numpoints=1, fancybox=True, loc='upper right')
    pylab.xlim(300, 1100)
    pylab.ylim(0, 4)
    ax.tick_params(axis='y', labelleft='off')
    ax.tick_params(axis='x', labelbottom='off') #labelsize='small')
    pylab.xticks(rotation=-45)
    #pylab.xlabel("Wavelength (nm)", fontsize='small')
    pylab.ylabel("F_lambda",  fontsize='small')    
    pylab.title("Example Stellar SEDs",  fontsize='small')

    # plot atmospheres
    ax = pylab.subplot(534)
    key = "Standard"
    ax.plot(atmos[key].wavelen, atmos[key].sb, 'b:', label=key)
    key = "X=1.0, H2O=0.7"
    ax.plot(atmos[key].wavelen, atmos[key].sb, 'g-', label=key)
    ax.tick_params(axis='y', labelsize='small')
    ax.tick_params(axis='x', labelbottom='off') #labelsize='small')
    pylab.xticks(rotation=-45)
    pylab.ylabel("Transmission",  fontsize='small')
    pylab.xlim(300, 1100)
    pylab.ylim(0, 1)
    #pylab.xlabel("Wavelength (nm)", fontsize='small')
    pylab.title("X=1.0",  fontsize='small')
    #pylab.legend(numpoints=1, fancybox=True, loc=(0.1, 0.05))
    ax = pylab.subplot(535)
    key = "Standard"
    ax.plot(atmos[key].wavelen, atmos[key].sb, 'b-', label=key)
    ax.tick_params(axis='y', labelleft='off')
    ax.tick_params(axis='x', labelbottom='off') # labelsize='small')
    pylab.xlim(300, 1100)
    pylab.ylim(0, 1)
    pylab.title("Standard atmosphere", fontsize='small')
    ax = pylab.subplot(536)
    key = "Standard"
    ax.plot(atmos[key].wavelen, atmos[key].sb, 'b:',  label=key)
    key = "X=1.5, H2O=1.4"
    ax.plot(atmos[key].wavelen, atmos[key].sb, 'g-', label=key)
    ax.tick_params(axis='y', labelleft='off') #labelsize='small')
    ax.tick_params(axis='x', labelbottom='off') #labelsize='small')
    pylab.xticks(rotation=-45)
    #pylab.xlabel("Wavelength (nm)", fontsize='small')        
    #pylab.ylabel("Transmission",  fontsize='small')
    pylab.xlim(300, 1100)
    pylab.title("X=1.5",  fontsize='small')
    #pylab.legend(numpoints=1, fancybox=True, loc=(0.1, 0.05))

    # plot system throughputs    
    ax = pylab.subplot(537)
    i = 0
    for f in filterlist:
        ax.plot(sys_std[f].wavelen, sys_std[f].sb, colors[i]+"-", label=f)
        ax.plot(sys_edge[f].wavelen, sys_edge[f].sb, colors[i]+":")
        i = color_counter_next(i)
    pylab.ylim(0, 1)
    ax.tick_params(axis='y', labelsize='small')
    ax.tick_params(axis='x', labelsize='small')
    pylab.xticks(rotation=-45)
    pylab.xlabel("Wavelength (nm)", fontsize='small')
    pylab.ylabel("Transmission",  fontsize='small')
    #pylab.legend(numpoints=1, fancybox=True, loc='upper right')
    pylab.title("Sys: Std + 1% shift",  fontsize='small')
    pylab.xlim(300, 1100)

    ax = pylab.subplot(538)
    i = 0
    for f in filterlist:
        ax.plot(sys_std[f].wavelen, sys_std[f].sb, colors[i]+"-", label=f)
        ax.plot(sys_edge[f].wavelen, sys_edge[f].sb, colors[i]+":")
        i = color_counter_next(i)
    pylab.ylim(0, 1)
    ax.tick_params(axis='y', labelleft='off') # labelsize='small')
    ax.tick_params(axis='x', labelsize='small')
    pylab.xticks(rotation=-45)
    pylab.xlabel("Wavelength (nm)", fontsize='small')
    #pylab.ylabel("Transmission",  fontsize='small')
    #pylab.legend(numpoints=1, fancybox=True, loc='upper right')
    pylab.title("Sys: Std + 1% shift",  fontsize='small')
    pylab.xlim(300, 1100)

    ax = pylab.subplot(539)
    i = 0
    for f in filterlist:
        ax.plot(sys_std[f].wavelen, sys_std[f].sb, colors[i]+"-", label=f)
        ax.plot(sys_edge[f].wavelen, sys_edge[f].sb, colors[i]+":")
        i = color_counter_next(i)
    pylab.ylim(0, 1)
    ax.tick_params(axis='y', labelleft='off') # labelsize='small')
    ax.tick_params(axis='x', labelsize='small')
    pylab.xticks(rotation=-45)
    pylab.xlabel("Wavelength (nm)", fontsize='small')
    #pylab.ylabel("Transmission",  fontsize='small')
    #pylab.legend(numpoints=1, fancybox=True, loc='upper right')
    pylab.title("Sys: Std + 1% shift",  fontsize='small')
    pylab.xlim(300, 1100)

    # Add the delta mags at the bottom.
    ax = pylab.subplot(5,3,13)
    i = 0
    xspace = [0, 1, 2, 3, 4, 5]
    for f in filterlist:
        key = 'X=1.0, std sys'
        ax.plot(xspace[i], mags_delta['blue'][key][f], colors[i]+"*")
        ax.plot(xspace[i], mags_delta['red'][key][f], colors[i]+"*")
        key = 'X=1.0, +1% sys'
        ax.plot(xspace[i], mags_delta['blue'][key][f], colors[i]+"o")
        ax.plot(xspace[i], mags_delta['red'][key][f], colors[i]+'o')
        i = color_counter_next(i)
    pylab.ylim(-0.05, 0.03)
    ax.tick_params(axis='y', labelsize='small')
    #pylab.yticks(ticks)
    ax.tick_params(axis='x', labelsize='small')
    pylab.xlim(-1, 6)
    pylab.grid(which='major')
    pylab.xticks(xspace, filterlist)
    pylab.title("Delta Mag (X=1.0)",  fontsize='small')
    
    ax = pylab.subplot(5,3,14)
    i = 0
    xspace = [0, 1, 2, 3, 4, 5]
    tmp = numpy.zeros(len(xspace))
    for f in filterlist:
        key = 'std atm, std sys'
        ax.plot(xspace[i],tmp[i], colors[i]+"*")
        ax.plot(xspace[i], tmp[i], colors[i]+"*")
        key = 'std atm, +1% sys'
        ax.plot(xspace[i], mags_delta['blue'][key][f], colors[i]+"o")
        ax.plot(xspace[i], mags_delta['red'][key][f], colors[i]+"o")
        i = color_counter_next(i)
    pylab.ylim(-0.05, 0.03)
    ax.tick_params(axis='y', labelleft='off')
    ax.tick_params(axis='x', labelsize='small')
    pylab.xlim(-1, 6)
    pylab.grid(which='major')
    pylab.xticks(xspace, filterlist)
    pylab.title("Delta Mag (Std)",  fontsize='small')
    
    ax = pylab.subplot(5,3,15)
    i = 0
    xspace = [0, 1, 2, 3, 4, 5]
    for f in filterlist:
        key = 'X=1.5, std sys'
        ax.plot(xspace[i], mags_delta['blue'][key][f], colors[i]+"*")
        ax.plot(xspace[i], mags_delta['red'][key][f], colors[i]+"*")
        key = 'X=1.5, +1% sys'
        ax.plot(xspace[i], mags_delta['blue'][key][f], colors[i]+"o")
        ax.plot(xspace[i], mags_delta['red'][key][f], colors[i]+'o')
        i = color_counter_next(i)
    pylab.ylim(-0.05, 0.03)
    pylab.xlim(-1, 6)
    pylab.grid(which='major')
    ax.tick_params(axis='x', labelsize='small')
    pylab.xticks(xspace, filterlist)
    ax.tick_params(axis='y', labelleft='off', labelright='on', labelsize='small')
    pylab.title("Delta Mag (X=1.5)",  fontsize='small')                       



if __name__ == "__main__":

    sys_std, sys_edge = read_hardware()
    atmos = read_atmos()
    total = combine_throughputs(atmos, sys_std, sys_edge)
    stars = read_seds(total)
    delta_mags = calc_mags(stars, total)
    make_figures(stars, sys_std, sys_edge, atmos, total, delta_mags)
    pylab.show()
    
