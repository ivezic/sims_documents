import numpy
import pylab
from SimpleSb import SimpleSb
from SedSets import SedSets
from GhostData import GhostData

filterlist = ('u', 'g', 'r', 'i', 'z', 'y4')
ss = SimpleSb(filterlist=filterlist)

#1nm 

ss.shift_throughputs(1.)
pylab.figure()
ss.plot_throughputs_diff()
pylab.title('1 nm shift')
pylab.savefig('through_shift1.png',format='png')

ghostfile = 'camera_ghosting_ff_calibbias_jdsu_lsstcone_121128'
filterlist = ('u', 'g', 'r', 'i', 'z','y4')
gd = GhostData(ghostData=ghostfile+'.fits')

#ok, let's do this at an r of 300, so index 60
pylab.figure()
i=60
for f in filterlist:
    pylab.plot(gd.direct[f][i].wavelen,gd.direct[f][i].phi-gd.flatghost[f][i].phi, label=f)

pylab.xlabel('Wavelength (nm)')
pylab.ylabel(r'$\Delta\Phi$ (nm$^{-1}$)')
pylab.title('Ghosting, r=%3.1f'%gd.radii[i])
pylab.legend(loc=4)
pylab.savefig('ghost_deltaphi.png', format='png')
