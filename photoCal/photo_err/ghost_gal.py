import numpy
import pylab
from GhostData import GhostData
from SedSets import SedSets

ghostfile = 'camera_ghosting_ff_calibbias_jdsu_lsstcone_121128'
filterlist = ('u', 'g', 'r', 'i', 'z','y4')
gd = GhostData(ghostData=ghostfile+'.fits')
graymags = {}
for f in filterlist:
    graymags[f] = gd.calc_grayscale(f)
    pylab.figure()
    pylab.plot(gd.radii, graymags[f])
    pylab.xlabel('Radius')
    pylab.ylabel('Delta Mag (mmag)')
    pylab.title('Grayscale Illumination Correction\n %s -- %s' %(ghostfile, f))
    pylab.savefig('%s_%s_dmag_ic.png' %(ghostfile, f), format='png')

sn = SedSets()
sn.read_eline_galaxy(redshifts=numpy.arange(0,2,.1))

directmags = {}
ghostmags = {}
dmags = {}

for f in filterlist:
    directmags[f] = numpy.zeros((len(sn.gallist), len(gd.radii)), 'float')
    ghostmags[f] = numpy.zeros((len(sn.gallist), len(gd.radii)), 'float')
    dmags[f] = numpy.zeros((len(sn.gallist), len(gd.radii)), 'float')
    for i,s in enumerate(sn.gallist):
        directmags[f][i], ghostmags[f][i], dmags[f][i] = gd.calc_mags(sn.gals[s], f)

# Set up color info for sn.
gi = directmags['g'][:,0] - directmags['i'][:,0]
ug = directmags['u'][:,0] - directmags['g'][:,0]

for f in filterlist:
    print 'Working on filter %s' %(f)
    # Write the (same as previously horrible) file
    outfilename = 'jdsu_%s_sn' %(f)                                                              
    outfile = open(outfilename, 'w')                                                                
    print >> outfile, '#', ghostfile                                                                
    print >>outfile, '# snname epoch z u-g g-r [primary mags] [flatfield/ghosting mags] [dmag]'                                                                                      
    #print >>outfile, '# radii:', gd.radii
    #for i,s in enumerate(sn.gallist):
    #    print >>outfile, s, sn.epochs[i], sn.redshifts[i], ug[i], gi[i], \
    #        directmags[f][i], ghostmags[f][i], dmags[f][i]

zs = sn.redshifts
idx = numpy.argsort(zs)
fig = pylab.figure()
sp = 0
for f in filterlist:    
    pylab.subplot(231+sp)
    sp+=1
    for r in [0, ]:
        pylab.plot(zs[idx], dmags[f][:,r][idx], marker='.', linestyle='-', label='r=%.0f' %(gd.radii[r]))
    #for r in [20, 30, 40, 50, 60]:
    #    pylab.plot(zs[idx], dmags[f][:,r][idx], marker='.', linestyle='', label='r=%.0f' %(gd.radii[r]))
    for r in [len(gd.radii)-1, ]:
        pylab.plot(zs[idx], dmags[f][:,r][idx], marker='.', linestyle='-', label='r=%.0f' %(gd.radii[r]))
    pylab.xlabel('z (redshift)')
    pylab.ylabel('Delta Mag (mmag)')
    if sp == 3:  pylab.legend( loc=3,numpoints=1, fancybox=True)
    pylab.title('%s' %(f))
    #pylab.savefig('gal_%s_%s_dmag_radius.png' %(ghostfile, f), format='png')

fig.tight_layout()
pylab.savefig('gal_ghosting_dmag.png',format='png')
#pylab.show()
