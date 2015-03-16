import os
import numpy
import pylab
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass
from lsst.sims.catalogs.measures.photometry.Sed import Sed
import plot_dmags as pd

def read_atmospheres(atmStrip, rootDir, atmtype):
    tmp = os.listdir(rootDir)
    atmfiles = []
    atmKeys = []
    for t in tmp:
        if t.startswith('am1.2') and t.endswith('.plt'):
            atmfiles.append(t)
            if (atmtype == 'aerosol') | (atmtype=='ozone'):
                atmKeys.append(float(t.split('_')[2].strip(atmStrip).strip('.plt')))
            else:
                atmKeys.append(float(t.split('_')[1].strip(atmStrip).strip('.plt')))
    atmKeys = numpy.array(atmKeys)
    atmDict = {}
    for k, filename in zip(atmKeys, atmfiles):
        atmDict[k] = Bandpass()
        atmDict[k].readThroughput(os.path.join(rootDir,filename))
    return atmDict, atmKeys


if __name__ == '__main__':
    # Read the water vapor atmospheres from disk
    pwv, pwv_values = read_atmospheres('pwv', '../water_sealevel', 'watervapor')
    # Read the aerosol atmospheres.
    aerosol, aerosol_values = read_atmospheres('avis', '../aerosol', 'aerosol')
    # Read the ozone atmospheres.
    ozone, ozone_values = read_atmospheres('oz', '../ozone', 'ozone')

    # Read the telescope hardware.
    filterlist = ('u', 'g', 'r', 'i', 'z', 'y3', 'y4')
    hardware = pd.read_hardware(shift_perc=None, flist=filterlist)

    # read the kurucz stars from disk
    stars, starlist, temperature, met, logg = pd.read_kurucz()
    

    atmosphere = [pwv, aerosol, ozone]
    atmosphereKeys = [pwv_values, aerosol_values, ozone_values]
    atmosphereName = ['Water Vapor (mm)', 'Aerosol (m - visibility)', 'Ozone (Dobson)']
    for atm, atmkeys, atmname in zip(atmosphere, atmosphereKeys, atmosphereName):
        # Set up a reference atmosphere/hardware .. use 'middle' value
        if ((len(atmkeys)%2) == 0):
            refkey = numpy.median(atmkeys[0:len(atmkeys)-1])
        else:
            refkey = numpy.median(atmkeys)
        nummags = len(starlist)
        iimags = {}
        for f in filterlist:
            iimags[f] = numpy.zeros([len(atmkeys), nummags], float)
        print numpy.shape(iimags['g'])
        # Calculate all the magnitudes.
        print atmkeys
        for i in range(len(atmkeys)):
            atmval = atmkeys[i]
            print atmval
            lsst = pd.combine_throughputs(atm[atmval], hardware)
            mags = pd.calc_adu(stars, starlist, lsst, filterlist)
            for f in filterlist:
                iimags[f][i] = numpy.copy(mags[f])
        #  Pull out the 'reference' magnitudes and colors.
        refidx = numpy.where(atmkeys==refkey)[0][0]
        refmags = {}
        for f in filterlist:
            lsst = pd.combine_throughputs(atm[refkey], hardware)
            mags = pd.calc_mags(stars, starlist, lsst, filterlist)
            refmags[f] = mags[f]
        refgi = refmags['g'] - refmags['i']
        gisort = numpy.argsort(refgi)
        refgi = numpy.copy(refgi[gisort])
        refmags = {}
        for f in filterlist:
            refmags[f] = numpy.copy(iimags[f][refidx])
        # subtract refmags from each line in iimags and convert to mmags
        for f in filterlist:
            for i in range(len(atmkeys)):
                iimags[f][i] = (iimags[f][i] - refmags[f])*1000.0
                iimags[f][i] = iimags[f][i][gisort]
        # Create figure showing changes in magnitude.
        atm_ext = numpy.array([(atmkeys.max() + (atmkeys.max() - atmkeys[len(atmkeys)-2])),], float)
        gatmkeys = numpy.concatenate((atmkeys, atm_ext))
        for f in filterlist:
            pylab.figure()
            vmax = iimags[f].max()
            vmin = iimags[f].min()
            if vmax < 1:
                vmax = 1.0
            if vmin > -1:
                vmin = -1.0
            pylab.pcolor(refgi, gatmkeys, iimags[f], vmin=vmin, vmax=vmax)
            cb = pylab.colorbar()
            cb.set_label(r'$\Delta$ %s (mmag)' %(f))
            pylab.xlim(-0.9, 1.5)
            pylab.xlabel('$g-i$')
            pylab.ylabel('%s' %(atmname))
            pylab.title('%s' %(atmname.split()[0]))
            figname = atmname.split()[0] + '_grid_' + f + '.png'
            pylab.savefig(figname, format='png')
        # And min/max values against median.
        metallicity = met[gisort]
        metcolors = ['c', 'c', 'b', 'g', 'y', 'r', 'm']
        metbinsize = abs(metallicity.min() - metallicity.max())/6.0
        metbins = numpy.arange(metallicity.min(), metallicity.max() + metbinsize, metbinsize)
        for f in filterlist:
            pylab.figure()
            for metidx in range(len(metbins)):
                condition =((metallicity>=metbins[metidx]) & (metallicity<=metbins[metidx]+metbinsize))
                mcolor = metcolors[metidx]
                pylab.plot(refgi[condition], iimags[f][0][condition], mcolor+'x')
                pylab.plot(refgi[condition], iimags[f][len(atmkeys)-1][condition], mcolor+'+')
            ymin, ymax = pylab.ylim()
            if ymin > -1:
                ymin = -1
            if ymax < 1:
                ymax = 1
            pylab.ylim(ymin, ymax)
            pylab.xlabel('$g-i$')
            pylab.ylabel(r"$\Delta$ %s (mmag)" %(f))
            pylab.title('%s band effect of min/max %s' %(f, atmname.split()[0]))
            pylab.grid()
            figname = atmname.split()[0] + '_maxmin_' +f + '.png'
            pylab.savefig(figname, format='png')
        # and change in all mags for g-i=0.5 object.
        colors = {'u':'c', 'g':'g', 'r':'r', 'i':'m', 'z':'k', 'y3':'y', 'y4':'burlywood'}
        pylab.figure()
        for f in filterlist:
            newrefmags = numpy.copy(iimags[f][0])            
            condition = (numpy.abs(refgi - 0.5) < 0.01)
            dmags_ave = numpy.zeros(len(atmkeys), float)
            dmags_err = numpy.zeros([2,len(atmkeys)], float)
            for i in range(len(atmkeys)):
                iimags[f][i] = iimags[f][i] - newrefmags
                tempdmags = iimags[f][i][condition]
                print '0.5', f, atmkeys[i], len(tempdmags)
                dmags_ave[i] = tempdmags.mean()
                dmags_err[0][i] = tempdmags.min()
                dmags_err[1][i] = tempdmags.max()
            pylab.errorbar(atmkeys, dmags_ave, dmags_err, elinewidth=1, linewidth=2,
                           color=colors[f], linestyle='-', marker='o', label=f)
        if atmname.startswith('Ozone'):
            pylab.legend(loc='lower left', numpoints=1, fancybox=True, shadow=True)
        else:
            pylab.legend(loc='upper left', numpoints=1, fancybox=True, shadow=True)
        atmkeys_step = atmkeys[1] - atmkeys[0]
        pylab.xlim(atmkeys.min()-atmkeys_step, atmkeys.max()+atmkeys_step)
        pylab.xlabel('%s' %(atmname))
        pylab.ylabel(r'$\Delta$ mag for g-i=0.5 (mmag)')
        pylab.title('%s' %(atmname.split()[0]))
        pylab.grid()
        figname = atmname.split()[0] + '_gradient1.png'
        pylab.savefig(figname, format='png')
        # and change in all mags for g-i=0.5 object
        pylab.figure()
        for f in filterlist:
            condition = (numpy.abs(refgi - (-0.5)) < 0.01)
            dmags_ave = numpy.zeros(len(atmkeys), float)
            dmags_err = numpy.zeros([2,len(atmkeys)], float)
            for i in range(len(atmkeys)):
                tempdmags = iimags[f][i][condition]
                print '-0.5', f, atmkeys[i], len(tempdmags)
                dmags_ave[i] = tempdmags.mean()
                dmags_err[0][i] = tempdmags.min()
                dmags_err[1][i] = tempdmags.max()
            pylab.errorbar(atmkeys, dmags_ave, dmags_err, elinewidth=1,  linewidth=2,
                           color=colors[f], linestyle='-', marker='o', label=f)
        if atmname.startswith('Ozone'):
            pylab.legend(loc='lower left', numpoints=1, fancybox=True, shadow=True)
        else:
            pylab.legend(loc='upper left', numpoints=1, fancybox=True, shadow=True)
        atmkeys_step = atmkeys[1] - atmkeys[0]
        pylab.xlim(atmkeys.min()-atmkeys_step, atmkeys.max()+atmkeys_step)
        pylab.xlabel('%s' %(atmname))
        pylab.ylabel(r'$\Delta$ mag for g-i=-0.5 (mmag)')
        pylab.title('%s' %(atmname.split()[0]))
        pylab.grid()
        figname = atmname.split()[0] + '_gradient2.png'
        pylab.savefig(figname, format='png')
