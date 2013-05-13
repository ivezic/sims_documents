import os
import numpy
import pylab


deg2rad = numpy.pi/180.0
rad2deg = 180.0/numpy.pi

colors = ('b', 'm', 'r', 'g', 'y', 'k')
def next_color(colorindex):
    colorindex = colorindex + 1
    if colorindex > len(colors):
        colorindex = 0
    return colorindex


def readModtranFile(filename):
    # read modtran files from David Burke
    # modtran file names are like Pachon_MODTRAN.12.7sc  .. 12=X=1.2
    # files have 2 header lines which include unknown info
    # then 2 lines of info on columns
    #  WAVELENGTH COMBIN    H2O   UMIX     O3  TRACE     N2    H2O MOLEC AER+CLD  HNO3 AER+CLD
    #   (NM)  TRANS  TRANS  TRANS  TRANS  TRANS   CONT   CONT   SCAT  TRANS  TRANS abTRNS
    # then continue until -9999 flag for wavelength at end of file
    # 
    # David Burke's code (MODTRAN_py4lj.txt) indicates: 
    # Read MODTRAN transmission templates formated by Modtran_ATMO_format.py.
    # Translation from MODTRAN to ATMO_fit:
    #   H2O = H2O TRANS * H2O SCAT  (but looks like H2O (CONT)?)
    #   O2  = UMIX * TRACE  = molecular absorption
    #   O3  = O3 = ozone
    #   Rayleigh = MOLEC SCAT = molecular scattering
    #   aerosol = Built-in power-law model. (but AER+CLD trans* abTrans it looks like)
    #        (although in these sims, =1.0 all the way through)

    # comb = combined version (i.e. absorption at that airmass for these components
    # but really, want to split out components at each airmass and fitting would use these
    # individual components to build an atmosphere.
    
    fin = open(filename, 'r')   
    rwavelen = []
    atm_comp = ('comb', 'H2O', 'O2', 'O3', 'rayleigh', 'aerosol')
    rtrans = {}
    for comp in atm_comp:
        rtrans[comp] = []
    i = 0
    for lines in fin:
        if i<4:
            if i == 1:
                values = lines.split()
                secz = 1/numpy.cos((float(values[2]))*deg2rad)
                secz = round(secz*10)/10.0
            i = i+1
            continue # skip over header lines
        values = lines.split()
        if (float(values[0]) < 0):
            break  # end of file marker (-9999.)
        rwavelen.append(values[0])
        rtrans['comb'].append(values[1])
        rtrans['H2O'].append(float(values[2])*float(values[7]))
        rtrans['O2'].append(float(values[3])*float(values[5]))
        rtrans['O3'].append(values[4])
        rtrans['rayleigh'].append(values[8])
        rtrans['aerosol'].append(float(values[9])*float(values[11]))
    fin.close()
    
    rwavelen = numpy.array(rwavelen, dtype=float)
    # wavelength in nanometers
    wavelen = numpy.arange(300.0, 1201.0, 1.0, dtype=float)
    trans = {}
    templates = {}
    for comp in atm_comp:
        rtrans[comp] = numpy.array(rtrans[comp], dtype=float)
        # resample transmission onto the same wavelength grid
        trans[comp] = numpy.interp(wavelen, rwavelen, rtrans[comp], left=0.0, right=0.0)
        templates[comp] = 1.0 - trans[comp]
    # fix aerosols to somewhat realistic values
    trans['aerosol'] = numpy.exp(-(0.039 * secz)*(wavelen/675.0)**-1.7)
    templates['aerosol'] = 1.0 - trans['aerosol']
    return secz, wavelen, trans, templates, atm_comp


def plotAtmos(secz, wavelen, trans, atm_comp, xlim=(300, 1100), newfig=True, savefig=False, figroot='atmos'):
    if newfig:
        pylab.figure()
    colorindex = 0
    transtotal = numpy.ones(len(wavelen), dtype='float')
    for comp in atm_comp:
        transtotal = transtotal * trans[comp]
        pylab.plot(wavelen, trans[comp], colors[colorindex], label=comp)
        colorindex = next_color(colorindex)
    pylab.plot(wavelen, transtotal, 'k:')
    pylab.legend(loc='lower right', numpoints=1, fancybox=True)
    pylab.xlim(xlim[0], xlim[1])
    pylab.xlabel("Wavelength (nm)")
    pylab.ylabel("Transmission")
    pylab.title("Airmass %.2f" %(secz))
    return

def plotTemplates(secz, wavelen, templates, atm_comp, xlim=(300, 1100), newfig=True, savefig=False, figroot='atmos'):
    if newfig:
        pylab.figure()
    colorindex = 0
    for comp in atm_comp:
        pylab.plot(wavelen, templates[comp], colors[colorindex])
        colorindex = next_color(colorindex)
    xloc = 0.8
    yloc = 0.4
    colorindex = 0        
    for comp in atm_comp:
        pylab.figtext(xloc, yloc, comp, color=colors[colorindex])
        colorindex = next_color(colorindex)
        yloc = yloc - 0.03
    pylab.xlim(xlim[0], xlim[1])
    pylab.xlabel("Wavelength (nm)")
    pylab.ylabel("Template")
    pylab.title("Airmass %.2f" %(secz))
    return



def plotAtmosRatio(seczlist, wavelen, atm_trans, atm_comp, xlim=(300, 1100),
                   newfig=True, savefig=False, figroot='atmos_ratio'):
    trans_ratio = numpy.zeros(len(wavelen), dtype='float')
    zref = seczlist[0]
    eps = 1e-30
    for comp in atm_comp: 
        tmp = numpy.where(atm_trans[zref]==0, eps, atm_trans[zref])
        pylab.figure()        
        for z in seczlist:
            trans_ratio = atm_trans[z][comp] / atm_trans[zref][comp]
            #trans_ratio = numpy.where(numpy.isnan(trans_ratio), 0.0, trans_ratio)
            #trans_ratio = numpy.where(numpy.isneginf(trans_ratio), 0.0, trans_ratio)
            #trans_ratio = numpy.where(numpy.isinf(trans_ratio), 0.0, trans_ratio)
            pylab.plot(wavelen, trans_ratio, label='X=%.2f' %(z))
        pylab.legend(loc='lower right', fancybox=True, numpoints=1)
        pylab.xlim(xlim[0], xlim[1])
        #pylab.ylim(0, 1)
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Ratio")
        pylab.title("Changes in transmission with airmass due to component %s" %(comp))
    return

def plotAtm(seczlist, wavelen, abs_atm, atm_comp, xlim=(300, 1100),
            newfig=True, savefig=False, figroot='atmos_piece'):
    for comp in atm_comp: 
        pylab.figure()
        for z in seczlist:
            pylab.plot(wavelen, abs_atm[z][comp], label='X=%.2f' %(z))
        #pylab.plot(wavelen, abs_atm[z][comp]/numpy.sqrt(2.5), label="2.5X * 1.3")
        pylab.legend(loc='lower right', fancybox=True, numpoints=1)
        pylab.xlim(xlim[0], xlim[1])
        #pylab.xlim(600, 1100)
        pylab.ylim(0, 1)
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Transmission")
        pylab.title("Transmission of component %s" %(comp))
    return    
        
def plotAbs(seczlist, wavelen, abs_atm, atm_comp, xlim=(300, 1100),
            newfig=True, savefig=False, figroot='atmos_template'):
    for comp in atm_comp: 
        pylab.figure()
        for z in seczlist:
            pylab.plot(wavelen, abs_atm[z][comp], label='X=%.2f' %(z))
        #pylab.plot(wavelen, abs_atm[z][comp]/numpy.sqrt(2.5), label="2.5X * 1.3")
        pylab.legend(loc='upper left', fancybox=True, numpoints=1)
        #pylab.xlim(xlim[0], xlim[1])
        pylab.xlim(600, 1100)
        pylab.ylim(0, 0.62)
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Absorption")
        pylab.title("Absorption due to component %s" %(comp))
    return
        

def buildAtmos(secz, wavelen, atmo_templates, atmo_ind, seczlist, C, xlim=[300, 1100]):
    # Burke paper says atmosphere put together as 
    # Trans_total (alt/az/time) = Tgray * (e^-Z*tau_aerosol(alt/az/t)) * 
    #         * (1 - C_mol * BP(t)/BPo * A_mol(Z))  -- line 2
    #         * (1 - sqrt(C_mol * BP(t)/BPo) * A_mol(Z))  -- 3
    #         * (1 - C_O3 * A_O3(A) )
    #         * (1 - C_H2O(alt/az/time) * A_H2O(Z))
    # Tau_aerosol = trans['aerosol'] ... but replace with power law (because here all 1's)
    #  typical power law index is about tau ~ lambda^-1
    # A_mol = trans['O2']
    
    # should interpolate between airmass steps of modtran output

    # secz = secz of this observation
    # wavelen / atmo_templates == building blocks of atmosphere, with seczlist / atmo_ind keys
    # C = coeffsdictionary = to, t1, t2, alpha0 (for aerosol), C_mol, BP, C_O3, C_H2O  values
    

    BP0 = 782 # mb
    trans_total = numpy.ones(len(wavelen), dtype='float')
    trans_total = trans_total * (1.0 - C['mol'] * C['BP']/BP0 * atmo_templates[secz]['rayleigh'])  \
        * ( 1 - numpy.sqrt(C['mol'] * C['BP']/BP0) * atmo_templates[secz]['O2']) \
        * ( 1 - C['O3'] * atmo_templates[secz]['O3']) \
        * ( 1 - C['H2O'] * atmo_templates[secz]['H2O'])

    aerosol = numpy.exp(-secz * (C['to'] + C['t1']*0.0 + C['t2']*0.0) * (wavelen/675.0)**C['alpha'])
    trans_total = trans_total * aerosol

    #trans_total = trans_total * (1 - atmo_templates[secz]['aerosol'])

    pylab.figure()
    pylab.subplot(212)
    colorindex = 0
    for comp in atmo_ind:
        if comp == 'aerosol':
            pylab.plot(wavelen, (1-aerosol), colors[colorindex], label='aerosol')
        else:
            pylab.plot(wavelen, atmo_templates[secz][comp], colors[colorindex], label='%s' %(comp))
        colorindex = next_color(colorindex)
    leg =pylab.legend(loc=(0.88, 0.3), fancybox=True, numpoints=1)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    pylab.xlim(xlim[0], xlim[1])
    pylab.ylim(0, 1.0)
    pylab.xlabel("Wavelength (nm)")
    pylab.subplot(211)
    pylab.plot(wavelen, atmo_trans[seczlist[0]]['comb'], 'r-', label='Standard')
    pylab.plot(wavelen, trans_total, 'k-', label='Observed')
    leg = pylab.legend(loc=(0.88, 0.2), fancybox=True, numpoints=1)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    pylab.xlim(xlim[0], xlim[1])
    pylab.ylim(0, 1.0)
    pylab.title("Atmosphere at X=%.2f" %(secz))

    return


def writeUsefulAtmos(secz, wavelen, trans):
    atmo = lm.teleThruput('NULL')
    atmo.wavelen = wavelen
    atmo.thruput = trans['comb']
    filename = 'atmos_' + '%.0f' %(secz*10.) + '.dat'
    atmo.writeThruput(filename)
    return


#######

if __name__ == '__main__' : 
    
    # these are the files that David Burke sent me for Pachon atmosphere
    filenames = os.listdir(".")
    atmofiles = []
    for f in filenames:
        if (f.startswith('Pachon_MODTRAN')) & (f.endswith('.7sc')):
            atmofiles.append(f)

    #print atmofiles

    wavelen = {}     
    atmo_trans = {}
    atmo_templates = {}
    seczlist = []
    for file in atmofiles:
        secz, wavelen[secz], atmo_trans[secz], atmo_templates[secz], atm_comp = readModtranFile(file)
        seczlist.append(secz)

    print seczlist

    atm_comp = ('comb', 'H2O', 'O2', 'O3', 'rayleigh', 'aerosol')
    atmo_ind = ('H2O', 'O2', 'O3', 'rayleigh', 'aerosol')
    
    # now have dictionaries of each component, with dictionary keys = secz (their airmass)
        
    # have a look at the plots?
    #for secz in seczlist:
    #    plotAtmos(secz, wavelen[secz], atmo[secz], atm_comp, newfig=True, savefig=False, figroot='atmos')    
    #    writeUsefulAtmos(secz, wavelen[secz], atmo[secz])

    for secz in seczlist[0], seczlist[10]:
        #plotAtmos(secz, wavelen[secz], atmo_trans[secz], atmo_ind, newfig=True, savefig=False, figroot='atmos')
        #plotTemplates(secz, wavelen[secz], atmo_templates[secz], atmo_ind, newfig=True, savefig=False,figroot='atmos')
        pass

    this_seczlist = [seczlist[0], seczlist[10], seczlist[15]]
    this_comp = ('rayleigh',)
    #plotAbs(this_seczlist, wavelen[this_seczlist[0]], atmo_templates, this_comp)
    #plotAtm(this_seczlist, wavelen[this_seczlist[0]], atmo_trans, this_comp)
    #plotAtmosRatio(this_seczlist, wavelen[this_seczlist[0]], atmo_trans, this_comp)

    z = seczlist[10]
    print secz
    C = {'O3':0.9, 'to':3.9/100.0, 't1':0.02/100.0, 't2':-0.03/100.0, 'alpha':-1.7, 
         'mol':0.96, 'BP':782, 'H2O':0.9}
    buildAtmos(z, wavelen[seczlist[0]], atmo_templates, atmo_ind, seczlist, C)
    pylab.show()
