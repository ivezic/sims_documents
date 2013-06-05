import numpy
import pylab
from FilterShift import FilterShift
from SedSets import SedSets

writeFiles = False
max_jitter = 1.65

# Set up throughputs package (to get default directory for baseline throughput curves for FilterShift). 
filterlist = ('u', 'g', 'r', 'i', 'z','y4')
fs = FilterShift(filterlist=filterlist, max_jitter=max_jitter)

# Read in Kurucz models.
kurucz = SedSets()
kurucz.read_kurucz()

# Set up dictionaries to store magnitudes. 
mags = {}
dmags = {}

# Calculate magnitudes for stars. 
for f in filterlist:
    mags[f] = numpy.zeros((len(kurucz.starlist), len(fs.radii)), 'float')
    dmags[f] = numpy.zeros((len(kurucz.starlist), len(fs.radii)), 'float')
    for i,s in enumerate(kurucz.starlist):
        mags[f][i], dmags[f][i] = fs.calc_mags(kurucz.stars[s], f)

# Pull out central radius / central bandpass information as 'baseline'. 
basemags = {}
for f in filterlist:
    basemags[f] = mags[f][:,0]

# Set up color info for stars.
gi = basemags['g'] - basemags['i']
ug = basemags['u'] - basemags['g']

if writeFiles:
    # Write out data. 
    for f in filterlist:
        print 'Working on filter %s' %(f)
        # Write the (same as previously horrible) file
        outfilename = 'filterShift_%s_stars' %(f)                                                              
        outfile = open(outfilename, 'w')                                                                
        print >> outfile, '# Filter Shifts'                                                                
        print >>outfile, '# starname temperature metallicity logg u-g g-r [mags] [dmags]'
        print >>outfile, '# radii:', fs.radii
        for i,s in enumerate(kurucz.starlist):
            print >>outfile, s, kurucz.temperature[i], kurucz.met[i], kurucz.logg[i], ug[i], gi[i], \
                mags[f][i], dmags[f][i]

# Make some plots. 
idx = numpy.argsort(gi)

# Plot shift as a function of color, for select radii.
for f in filterlist:
    pylab.figure()
    for r in [0, ]:
        pylab.plot(gi[idx], dmags[f][:,r][idx], marker='.', linestyle='-', label='r=%.0f' %(fs.radii[r]))
    mid = len(fs.radii)/2    
    for r in [mid-20, mid, mid+20]:
        pylab.plot(gi[idx], dmags[f][:,r][idx], marker='.', linestyle='', label='r=%.0f' %(fs.radii[r]))
    for r in [len(fs.radii)-1, ]:
        pylab.plot(gi[idx], dmags[f][:,r][idx], marker='.', linestyle='-', label='r=%.0f' %(fs.radii[r]))
    pylab.xlabel('g-i')
    pylab.ylabel('Delta Mag (mmag)')
    pylab.legend(loc=(0.93, 0.2), fontsize='smaller', numpoints=1, fancybox=True)
    pylab.title('%s -- %s' %('filter shift + jitter', f))
    pylab.savefig('%s_%s_dmag_color.png' %('filtershift', f), format='png')

# Plot shift, for a particular g-i value, as a function of radius.
condition1 = (numpy.abs(gi-0.5)<0.02)
condition2 = (numpy.abs(gi+0.5)<0.02)
for f in filterlist:
    pylab.figure()
    pylab.plot(fs.radii, numpy.mean(numpy.abs(dmags[f][condition1]), 0), color='g', marker='o', label='g-i~0.5')
    pylab.plot(fs.radii, numpy.mean(numpy.abs(dmags[f][condition2]), 0), color='r', marker='o', label='g-i~-0.5')
    pylab.xlabel('Radius (mm)')
    pylab.ylabel('Delta Mag (mmag)')
    pylab.title('%s -- %s' %('filter shift + jitter', f))
    pylab.legend(loc=(0.93, 0.2), fontsize='smaller', numpoints=1, fancybox=True)
    pylab.savefig('%s_%s_dmag_radius.png' %('filtershift', f), format='png')

pylab.show()

