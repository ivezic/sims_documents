import numpy
import pylab
from lsst.sims.atmosphere.clouds.Arma.ArmaSf import ArmaSf
from lsst.sims.atmosphere.clouds.Arma.Clouds import Clouds


def plotsq(x0,y0, len, linewidth=3):
    pylab.plot([x0,x0+len],[y0,y0],'k',linewidth=linewidth)
    pylab.plot([x0,x0+len],[y0+len,y0+len],'k',linewidth=linewidth)
    pylab.plot([x0,x0],[y0,y0+len],'k',linewidth=linewidth)
    pylab.plot([x0+len,x0+len],[y0,y0+len],'k',linewidth=linewidth)
    return


cloud = Clouds()
armasf = ArmaSf()

# Set up lambda_p / lambda_avg / lambda_s (lambda_s = anything, but just much smaller than lambda_*)
lambda_p = 500.
lambda_avg = 300.
lambda_s = 2.

# Set kappa and c values
num_clouds = 1
c = 0.55
kappa = 0.5
#kappa = numpy.random.normal(loc=0.5, scale=0.3, size=num_clouds)
#kappa = numpy.where(kappa<0, 0, kappa)

fov=1.8*2
pixscale=40.

sftheta, sfsf = armasf.CloudSf(lambda_p, lambda_avg, lambda_s, kappa, c)
cloud.makeCloudImage(sftheta, sfsf, kappa, fov=fov, pixscale=pixscale, oversample=1.0, verbose=True)

pylab.imshow(cloud.cloudimage, extent=(0,fov,0,fov))
cb = pylab.colorbar()
pylab.xlabel('degrees')
pylab.ylabel('degrees')
pylab.title('Simulated Cloud Image')




plotsq(1,2,3.6/5.)
plotsq(1,2,3.6/5./3)

pylab.xlim([0,fov])
pylab.ylim([0,fov])

cb.set_label('extinction (mags)')
pylab.savefig('cloud.eps')
