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
resolution=1.

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
ss.plot_throughputs_diff()
#file = open('r_phi.dat','w')
#print >>file, '#wavelength  phi'
#for i in numpy.arange(numpy.size(ss.basebp['r'].wavelen)): print >>file, ss.basebp['r'].wavelen[i], ss.basebp['r'].phi[i]
#file.close()

pylab.savefig('gal_delta_through.png',format='png')

# Read in Kurucz model stars. 
gal = SedSets()
gal.read_eline_galaxy(redshifts=numpy.arange(0,2.1,.05))

# Set up dictionaries to store magnitudes. 
mags_base = {}
mags_new = {}
dmags = {}

# Calculate magnitudes. 
for f in filterlist:
    mags_base[f] = numpy.zeros(len(gal.gallist), 'float')
    mags_new[f] = numpy.zeros(len(gal.gallist), 'float')
    dmags[f] = numpy.zeros(len(gal.gallist), 'float')
    for i,s in enumerate(gal.gallist):
        mags_base[f][i], mags_new[f][i], dmags[f][i] = ss.calc_mags(gal.gals[s], f)
    

gi = mags_base['g'] - mags_base['i']
zs = gal.redshifts

# Plot changes in magnitude (due to shifting / downgrading resolution of phi)
fig = pylab.figure()
i=0
for f in filterlist:
    pylab.subplot(321+i)
    i+=1
    pylab.plot(zs, dmags[f], marker='.', linestyle='')
    pylab.xlabel('z (redshift)')
    pylab.ylabel('Delta Mag (mmag)')
    pylab.axhline(y=10, color='r', linestyle='-')
    pylab.axhline(y=-10, color='r', linestyle='-')
    limits = pylab.axis()
    pylab.ylim([numpy.min([limits[2],-13]),numpy.max([limits[3],13])])
    figtitle = '%s' %(f)
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
    
#pylab.show()

fig.tight_layout()
pylab.savefig('simple_gal_%1.0f.png'%resolution, format='png')


pylab.figure()
gal.gals['0.0'].flambdaTofnu()
pylab.plot(gal.gals['0.0'].wavelen, gal.gals['0.0'].fnu)
pylab.ylabel(r'f$_\nu$')
pylab.xlabel(r'wavelength (nm)')
pylab.title('Starburst Galaxy')
pylab.xlim([100,1100])
#pylab.ylim([0,30])
pylab.savefig('galaxy_sed.png',format='png')


#plot u-band
pylab.figure()
gal.gals['1.6'].flambdaTofnu()
pylab.plot(gal.gals['1.6'].wavelen, gal.gals['1.6'].fnu/1e15, label='Galaxy')
pylab.ylabel(r'f$_\nu$, $\Delta\Phi$ (nm$^{-1}$)')
pylab.xlabel(r'wavelength (nm)')
pylab.xlim([300,450])
pylab.ylim([-.002,0.004])
pylab.plot(ss.basebp['u'].wavelen, ss.basebp['u'].phi-ss.newbp['u'].phi, label=r'$\Delta\Phi$')
pylab.legend()
pylab.savefig('gal_phi_z1.6.png',format='png')

exit()
    

