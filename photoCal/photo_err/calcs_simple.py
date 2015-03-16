import numpy
import pylab
from SimpleSb import SimpleSb 
from SedSets import SedSets

# Set up SimpleSb (which can shift and de-resolve bandpasses).
# Read in baseline throughput curves. 
filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
ss = SimpleSb(filterlist=filterlist)

shift=True
degrade_resolution = False

#[optional]
# Add shift to bandpasses. 
# Shifts in nanometers - can use dictionary or single number. 
if shift:
    shift = 1. 
    #shift = {'u':2., 'g':3., 'r':3., 'i':4., 'z':4., 'y':5} 
    ss.shift_throughputs(shift)
    # Preliminary tests indicate that a shift of 1 nm already creates ~5 mmag of error in MS stars in rigzy filters, >20mmag in ug.

#[optional, although must do this or shift bandpass]
# Deresolve bandpasses
# 'resolution' is in nanometers, and can use dictionary or single number
if degrade_resolution:
    resolution = 5.
    ss.deres_throughputs(resolution)
    # Some preliminary tests indicate that resolution =5 nm creates < 1 mmag error for MS stars

# Let's plot the bandpasses, for fun. 
ss.plot_throughputs()

# Read in Kurucz model stars. 
kurucz = SedSets()
kurucz.read_kurucz()

# Set up dictionaries to store magnitudes. 
mags_base = {}
mags_new = {}
dmags = {}

# Calculate magnitudes. 
for f in filterlist:
    mags_base[f] = numpy.zeros(len(kurucz.starlist), 'float')
    mags_new[f] = numpy.zeros(len(kurucz.starlist), 'float')
    dmags[f] = numpy.zeros(len(kurucz.starlist), 'float')
    for i,s in enumerate(kurucz.starlist):
        mags_base[f][i], mags_new[f][i], dmags[f][i] = ss.calc_mags(kurucz.stars[s], f)
    
gi = mags_base['g'] - mags_base['i']

# Plot changes in magnitude (due to shifting / downgrading resolution of phi)
for f in filterlist:
    pylab.figure()
    pylab.plot(gi, dmags[f], marker='.', linestyle='')
    pylab.xlabel('g-i')
    pylab.ylabel('Delta Mag (mmag)')
    figtitle = 'Filter %s' %(f)
    if shift:
        if isinstance(shift, dict):
            figtitle += ': shifted %.1f nm' %(shift[f])
        else:
            figtitle += ': shifted %.1f nm' %(shift)
    if degrade_resolution:
        if isinstance(resolution, dict):
            figtitle += ': res %.1f nm' %(resolution[f])
        else:
            figtitle += ': res %.1f nm' %(resolution)
    pylab.title(figtitle)
    pylab.savefig('simple_%s.png' %(f), format='png')
pylab.show()

exit()
    

