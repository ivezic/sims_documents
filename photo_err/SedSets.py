"""
Reads and holds information about large groups of SEDs, as they are expected to reside in the accompanying data
directories.
"""

import os
import numpy
from copy import deepcopy
from lsst.sims.catalogs.measures.photometry.Sed import Sed

class SedSets:
    """ Set up a dictionary of a bunch of seds."""
    def __init__(self):
        return

    def read_kurucz(self, dataDir='data/kurucz_r', g40_only=True):
        """Read selected kurucz seds. Stores metallicity, logg and temperature information as well."""
        # read kurucz model MS, [g40_only?] stars SEDs
        allfilelist = os.listdir(dataDir)
        starlist = []
        # make preliminary cut for ms, g40 stars.
        if g40_only:
            for filename in allfilelist:
                if filename[-3:] == 'g40':
                    starlist.append(filename)
        else:
            starlist = allfilelist
        atemperature = []
        amet = []
        alogg = []
        starlist2 = []
        # Make secondary cut for stars with temperature > 4000 deg (kurucz models are unreliable below 4000K).
        for s in starlist:
            tmp = s.split('_')
            met = float(tmp[0][2:])
            if tmp[0][1] == 'm':
                met = -1 * met        
            met = met/10.0
            temperature = float(tmp[1][:5])
            logg = float(tmp[2][1:])
            logg = logg / 10.0
            if (temperature > 4000.0):
                amet.append(met)
                atemperature.append(temperature)
                alogg.append(logg)
                starlist2.append(s)
        # Update what we're keeping for star information. 
        self.temperature = numpy.array(atemperature)
        self.met = numpy.array(amet)
        self.logg = numpy.array(alogg)
        self.starlist = starlist2
        # And now, actually read the stars SEDS from disk.
        self.stars = {}
        for s in self.starlist:
            self.stars[s] = Sed()
            self.stars[s].readSED_flambda(os.path.join(dataDir,s))
        print "# Read %d MS stars from %s" %(len(self.starlist), dataDir)
        return


    def read_sn(self, dataDir='data/sn', epochs=['10', '15', '20'], redshifts=numpy.arange(0, 1.0, 0.1)):
        """Reads SN at a variety of epochs and then redshifts to the desired range of z. Stores
        seds as well as epochs and redshifts. """
        # Read sn spectra and redshift
        allfilelist = os.listdir(dataDir)
        snlist = []
        # Pull out the filenames we want (which match epochs and SNIa template)
        for filename in allfilelist:
            if filename.endswith('.dat') & filename.startswith('sn1a_'):
                #snlist.append(filename)
                tmp = filename.split('_')
                epoch = tmp[1].split('.')[0]
                if epoch in epochs:
                    snlist.append(filename)
        # Read base SEDs for these epochs.
        sns_base = {}
        for sn, epoch in zip(snlist, epochs):
            sns_base[epoch] = Sed()
            sns_base[epoch].readSED_flambda(os.path.join(dataDir, sn))
        # Then redshift to build stored set of SEDs.
        self.sns = {}
        self.snlist = []
        for e in epochs:
            for z in redshifts:
                sn_name = "%d_%.1f" %(int(e), z)
                wavelen, flambda = sns_base[e].redshiftSED(z, wavelen=sns_base[e].wavelen,
                                                           flambda=sns_base[e].flambda)
                self.sns[sn_name] = Sed(wavelen=wavelen, flambda=flambda)
                self.snlist.append(sn_name)
        print "# Generated %d sn's at redshifts between %f and %f for epochs %s" %(len(self.snlist), \
                                                                                   redshifts.min(), \
                                                                                   redshifts.max(), epochs)
        self.epochs = epochs
        self.redshifts = redshifts
        return


    def read_whitedwarf(self, dataDir='data/white_dwarfs_r'):
        # read white dwarf bergeron models
        # get the H dwarfs
        Hdir = os.path.join(dataDir, "H")
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
        Hedir = os.path.join(dataDir, "He")
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
        self.temperatures = numpy.array(temperatures)
        self.loggs = numpy.array(loggs)
        self.hlist = hlist
        self.helist = helist
        self.wdlist = hlist + helist
        self.wds = {}
        for w in self.wdlist:
            self.wds[w] = Sed()
            if w in hlist:
                self.wds[w].readSED_flambda(os.path.join(Hdir, w))
            if w in helist:
                self.wds[w].readSED_flambda(os.path.join(Hedir, w))
        return

