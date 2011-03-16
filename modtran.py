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
#    trans['aerosol'] = trans['aerosol'] + (0.98-trans['aerosol'][wavelen==800])
    templates['aerosol'] = 1.0 - trans[comp]
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



def plotAtmosRatio(secz_a, wavelen_a, trans_a, secz_b, wavelen_b, trans_b, atm_comp, xlim=(300, 1000), newfig=True, savefig=False, figroot='atmos_ratio'):
    zeroes = numpy.zeros((len(wavelen_a),), dtype=float)
    trans_ratio = {}
    colorindex= 0
    for comp in atm_comp: #atm_comp:
        trans_ratio[comp] = trans_b[comp] / trans_a[comp]
        trans_ratio[comp] = numpy.where(numpy.isnan(trans_ratio[comp]), zeroes, trans_ratio[comp])
        trans_ratio[comp] = numpy.where((numpy.isneginf(trans_ratio[comp]) | numpy.isinf(trans_ratio[comp])), zeroes, trans_ratio[comp])
        pylab.figure()
        pylab.plot(wavelen_a, trans_a[comp], 'r:')
        pylab.plot(wavelen_a, trans_b[comp], 'b:')
        pylab.plot(wavelen_a, trans_ratio[comp], 'k--')
        pylab.figtext(0.8, 0.4, comp)
        pylab.xlim(xlim[0], xlim[1])
        pylab.ylim(0, 1.2)
        pylab.xlabel("Wavelength (nm)")
        pylab.ylabel("Ratio of Transmission (%.2f) and Transmission (%.2f)" %(secz_b, secz_a))
    return
        


def buildAtmos(secz, wavelenDict, transDict, comp_ratios, atmo_ind, seczlist):
    # wavelenDic and transDict = dictionaries of all secz, containing atmosphere components
    # Burke paper says atmosphere put together as 
    # Trans_total (alt/az/time) = Tgray * (e^-Z*tau_aerosol(alt/az/t)) * 
    #         * (1 - C_mol * BP(t)/BPo * A_mol(Z))  -- line 2
    #         * (1 - sqrt(C_mol * BP(t)/BPo) * A_mol(Z))  -- 3
    #         * (1 - C_O3 * A_O3(A) )
    #         * (1 - C_H2O(alt/az/time) * A_H2O(Z))
    # Tau_aerosol = trans['aerosol'] ... but replace with power law (because here all 1's)
    #  typical power law index is about tau ~ lambda^-1
    #    total absorption should also be proportional to airmass (right .. Z factor .. ok)
    # A_mol = trans['O2']
    
    # but Burke presentation from Calibration Boston 2009 says
    #  C_mol = C_BP but A_mol = trans_Rayleigh (trans_ray?) in line 2
    #  C_mol = C_BP but A_mol = trans_mol in line 3 (with sqrt) 
    #     but where is trans_O2???
    # also, when he was checking trans_comb, just multiplied all of these components together...
    # trans_total = 
    # should interpolate between airmass steps of modtran output
    trans_total = numpy.ones((len(allwavelen[secz]),))
    trans_total = (numpy.exp(-secz*alltrans[secz]['aerosol'])) #######
    
    for comp in atmo_ind:
        trans_total = trans_total *  alltrans[secz][comp]
    pylab.figure()
    #pylab.plot(allwavelen[secz], trans_total, 'k:')
    #pylab.plot(allwavelen[secz], alltrans[secz]['comb'], 'r:')
    
    # should try building an atmosphere according to above formula
    trans_total = numpy.ones((len(allwavelen[secz]),))
    C_BP = 0.9
    C_O3 = 0.9
    C_H2O = 0.9

    # but for now, let's cheat. 
    # the old atmos.dat file comes very close to secz=1.2, seczlist[2] file

    # just divide out water for low water vapor, and then add it back in.
    atmo = lm.teleThruput('NULL')
    atmo.wavelen = allwavelen[secz]
    atmo.thruput = alltrans[secz]['comb']
    waterthruput = alltrans[secz]['H2O']
    atmo.thruput = atmo.thruput/waterthruput
    
    waterthruput = alltrans[seczlist[0]]['H2O']
    atmo.thruput = atmo.thruput*waterthruput
    atmo.writeThruput("Atmo_12_lowater.dat")
    pylab.plot(atmo.wavelen, atmo.thruput, 'k')
    pylab.figtext(0.2, 0.38, 'X=1.2', color='k')
    pylab.figtext(0.2, 0.34, '~-30%s H2O' %('%'), color='k')

    atmo.thruput = atmo.thruput/waterthruput
    waterthruput = alltrans[seczlist[6]]['H2O']
    atmo.thruput = atmo.thruput * waterthruput
    atmo.writeThruput("Atmo_12_hiwater.dat")
    pylab.plot(atmo.wavelen, atmo.thruput, 'g')
    pylab.figtext(0.2, 0.3, '~+30%s H2O' %('%'), color='g')
    
    y3 = lm.teleThruput("/Users/ljones/work/filters/thruputs/final_y3")
    y4 = lm.teleThruput("/Users/ljones/work/filters/thruputs/final_y4")
    pylab.plot(y3.wavelen, y3.thruput*2., 'm:')
    pylab.plot(y4.wavelen, y4.thruput*2., 'b:')

    pylab.xlabel("Wavelength (Angstrom)")
    pylab.ylabel("Thruput, Percent")
    pylab.xlim(800, 1100)
    pylab.ylim(0.4, 1)

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
        plotAtmos(secz, wavelen[secz], atmo_trans[secz], atmo_ind, newfig=True, savefig=False, figroot='atmos')
        #plotTemplates(secz, wavelen[secz], atmo_templates[secz], atmo_ind, newfig=True, savefig=False, figroot='atmos')

    #plotAtmosRatio(seczlist[0], wavelen[seczlist[0]], atmo_trans[seczlist[0]], seczlist[15], wavelen[seczlist[15]], atmo_trans[seczlist[10]], atmo_ind)

    secz = seczlist[2]
    print secz
    #buildAtmos(secz, wavelen, atmo, atm_comp, seczlist)

    pylab.show()
