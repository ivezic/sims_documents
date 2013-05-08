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
sn.read_sn()

directmags = {}
ghostmags = {}
dmags = {}

for f in filterlist:
    directmags[f] = numpy.zeros((len(sn.snlist), len(gd.radii)), 'float')
    ghostmags[f] = numpy.zeros((len(sn.snlist), len(gd.radii)), 'float')
    dmags[f] = numpy.zeros((len(sn.snlist), len(gd.radii)), 'float')
    for i,s in enumerate(sn.snlist):
        directmags[f][i], ghostmags[f][i], dmags[f][i] = gd.calc_mags(sn.sns[s], f)

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
    #for i,s in enumerate(sn.snlist):
    #    print >>outfile, s, sn.epochs[i], sn.redshifts[i], ug[i], gi[i], \
    #        directmags[f][i], ghostmags[f][i], dmags[f][i]

idx = numpy.argsort(gi)

for f in filterlist:
    pylab.figure()
    for r in [0, ]:
        pylab.plot(gi[idx], dmags[f][:,r][idx], marker='.', linestyle='-', label='r=%.0f' %(gd.radii[r]))
    for r in [20, 30, 40, 50, 60]:
        pylab.plot(gi[idx], dmags[f][:,r][idx], marker='.', linestyle='', label='r=%.0f' %(gd.radii[r]))
    for r in [len(gd.radii)-1, ]:
        pylab.plot(gi[idx], dmags[f][:,r][idx], marker='.', linestyle='-', label='r=%.0f' %(gd.radii[r]))
    pylab.xlabel('g-i')
    pylab.ylabel('Delta Mag (mmag)')
    pylab.legend(loc=(0.93, 0.2),  numpoints=1, fancybox=True)
    pylab.title('%s -- %s' %(ghostfile, f))
    pylab.savefig('sn_%s_%s_dmag_radius.png' %(ghostfile, f), format='png')
#pylab.show()
