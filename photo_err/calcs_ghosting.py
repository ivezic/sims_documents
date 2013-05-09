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

kurucz = SedSets()
kurucz.read_kurucz()

directmags = {}
ghostmags = {}
dmags = {}

# Stars. 
for f in filterlist:
    directmags[f] = numpy.zeros((len(kurucz.starlist), len(gd.radii)), 'float')
    ghostmags[f] = numpy.zeros((len(kurucz.starlist), len(gd.radii)), 'float')
    dmags[f] = numpy.zeros((len(kurucz.starlist), len(gd.radii)), 'float')
    for i,s in enumerate(kurucz.starlist):
        directmags[f][i], ghostmags[f][i], dmags[f][i] = gd.calc_mags(kurucz.stars[s], f)

# Set up color info for stars.
gi = directmags['g'][:,0] - directmags['i'][:,0]
ug = directmags['u'][:,0] - directmags['g'][:,0]

for f in filterlist:
    print 'Working on filter %s' %(f)
    # Write the (same as previously horrible) file
    outfilename = 'jdsu_%s_stars' %(f)                                                              
    outfile = open(outfilename, 'w')                                                                
    print >> outfile, '#', ghostfile                                                                
    print >>outfile, '# starname temperature metallicity logg u-g g-r [primary mags] [flatfield/ghosting mags] [dmag]'                                                                                      
    print >>outfile, '# radii:', gd.radii
    for i,s in enumerate(kurucz.starlist):
        print >>outfile, s, kurucz.temperature[i], kurucz.met[i], kurucz.logg[i], ug[i], gi[i], \
            directmags[f][i], ghostmags[f][i], dmags[f][i]

# Plot dmags as a function of color, for specific radii.
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
    pylab.legend(loc=(0.93, 0.2), fontsize='smaller', numpoints=1, fancybox=True)
    pylab.title('%s -- %s' %(ghostfile, f))
    pylab.savefig('%s_%s_dmag_color.png' %(ghostfile, f), format='png')

# Plot shift, for a particular g-i value, as a function of radius.             
condition1 = (numpy.abs(gi-0.5)<0.02)
condition2 = (numpy.abs(gi+0.5)<0.02)
for f in filterlist:
    pylab.figure()
    pylab.plot(gd.radii, numpy.mean(dmags[f][condition1], 0), color='g', marker='o', label='g-i~0.5')
    pylab.plot(gd.radii, numpy.mean(dmags[f][condition2], 0), color='r', marker='o', label='g-i~-0.5')
    pylab.xlabel('Radius (mm)')
    pylab.ylabel('Delta Mag (mmag)')
    pylab.title('%s -- %s' %(ghostfile, f))
    pylab.legend(loc=(0.93, 0.2), fontsize='smaller', numpoints=1, fancybox=True)
    pylab.savefig('%s_%s_dmag_radius.png' %(ghostfile, f), format='png')

pylab.show()

