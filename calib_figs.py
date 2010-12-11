import numpy
import pylab

skysize = 250
psf = 10  # psf of a 'star' at center of fov
sensitivity_var = 0.002  # amount of variation around 1 (flat random) for pixel sensitivity
pixscale_change = -0.005 # fraction change pixel scale at center to outer edges. scale < 0 -> more flux at edges.
ghostscale = 1e-3  # fraction of total counts that are added in ghost pupil.
domegradient = 0.01  # fraction of change in total counts due to the dome screen not being illuminated evenly.

def make_pixgrid():
    # Define the pixel grid - this is the grid in the camera, where pixel distortion can occur,
    # adding sensitivity variations which are flat random, with mean of 1, +/- sensitivity_var.
    pixgrid = numpy.random.random([skysize, skysize]) * sensitivity_var * 2 + (1-sensitivity_var)
    pixgrid = pixgrid / pixgrid.mean()
    # hey, let's add a dust ring
    dustringx = 120
    dustringy = 80
    dustringrad_inner = 5
    dustringrad_outer = 9
    for i in range(dustringx-dustringrad_outer, dustringx+dustringrad_outer):
        for j in range(dustringy-dustringrad_outer, dustringy+dustringrad_outer):
            dist = numpy.sqrt((i-dustringx)**2+(j-dustringy)**2)
            if (dist < dustringrad_outer) & (dist > dustringrad_inner):
                if (i>0) & (i<skysize-1) & (j>0) & (j<skysize-1):
                    pixgrid[i][j] = pixgrid[i][j] - 2.5*sensitivity_var
                    if pixgrid[i][j] < 0:
                        pixgrid[i][j] = 0.001                    
    return pixgrid

def add_bkgd(pixgrid, bkgdval):
    pixgrid = pixgrid/pixgrid.mean()*bkgdval
    return pixgrid

def add_star(star, rawdetector):
    stargrid = numpy.zeros([skysize, skysize], 'float')
    # star = a list or tuple containing starx, stary, starflux
    starx = star[0]
    stary = star[1]
    starflux = star[2]
    # calculate star's extent at this X/Y by taking "PSF" = psf pix at center of fov
    pixcenX = skysize-1
    pixcenY = skysize-1 
    totaldistance = numpy.sqrt(pixcenX**2+pixcenY**2)
    distance = numpy.sqrt((starx-pixcenX)**2 + (stary-pixcenY)**2)
    starpsf = (1 + pixscale_change*(distance/totaldistance))*psf/2.0
    starflux_perpix = starflux / (starpsf**2*numpy.pi)
    starx0 = int(numpy.floor(starx - starpsf))
    starx1 = int(numpy.ceil(starx + starpsf))
    stary0 = int(numpy.floor(stary - starpsf))
    stary1 = int(numpy.ceil(stary + starpsf))
    #print starpsf, starx0, starx1, stary0, stary1
    starpix = 0
    for i in range(starx0, starx1):
        for j in range(stary0, stary1):
            # Send pixgrid as raw detector sensitivity at star's location (no bkdg added)!!
            if (numpy.sqrt((i-starx)**2+(j-stary)**2) <= starpsf):
                if (i>0) & (i<skysize-1) & (j>0) & (j<skysize-1):
                    stargrid[i][j] = stargrid[i][j] + rawdetector[i][j]*starflux
                    starpix = starpix + 1
    stargrid = stargrid / starpix
    return stargrid

def measure_star(pixgrid, star):
    starx = star[0]
    stary = star[1]
    starflux = star[2]
    # calculate star's extent at this x/y
    pixcenX = skysize-1
    pixcenY = skysize-1 
    totaldistance = numpy.sqrt(pixcenX**2+pixcenY**2)
    distance = numpy.sqrt((starx-pixcenX)**2 + (stary-pixcenY)**2)
    starpsf = (psf + psf * pixscale_change*distance/totaldistance)/2.0
    starx0 = int(numpy.floor(starx - starpsf))
    starx1 = int(numpy.ceil(starx + starpsf))
    stary0 = int(numpy.floor(stary - starpsf))
    stary1 = int(numpy.ceil(stary + starpsf))
    # measure pixel values within this aperture
    starmeas = 0
    starpix = 0
    for i in range(starx0, starx1):
        for j in range(stary0, stary1):
            if (numpy.sqrt((i-starx)**2+(j-stary)**2) <= starpsf):
                starmeas = starmeas + pixgrid[i][j]
                starpix = starpix + 1
    # get estimate for background here
    bkgdextent = 4
    bkgdpix = 0
    bkgd = 0
    for i in range(starx0 - bkgdextent, starx0-1):
        for j in range(stary0 - bkgdextent, stary0 + bkgdextent):
            if (i>0) & (i<skysize-1) & (j>0) & (j<skysize-1):
                bkgd = bkgd + pixgrid[i][j]
                bkgdpix = bkgdpix + 1
    for i in range(starx1+1, starx1+bkgdextent):
        for j in range(stary0 - bkgdextent, stary0 + bkgdextent):
            if (i>0) & (i<skysize-1) & (j>0) & (j<skysize-1):
                bkgd = bkgd + pixgrid[i][j]
                bkgdpix = bkgdpix + 1
    for i in range(starx0+1, starx1-1):
        for j in range(stary0-bkgdextent, stary0-1):
            if (i>0) & (i<skysize-1) & (j>0) & (j<skysize-1):
                bkgd = bkgd + pixgrid[i][j]
                bkgdpix = bkgdpix + 1
    for i in range(starx0+1, starx1-1):
        for j in range(stary1+1, stary1+bkgdextent):
            if (i>0) & (i<skysize-1) & (j>0) & (j<skysize-1):
                bkgd = bkgd + pixgrid[i][j]
                bkgdpix = bkgdpix + 1
    bkgd = bkgd / bkgdpix
    starmeas = starmeas - (bkgd*starpix)
    print starmeas, bkgd, starmeas + (bkgd*starpix)
    return starmeas, bkgd
        
def add_domewvariance(pixgrid, domeval):
    # Add a gradient across the field of view to represent variability in the dome screen illumination.
    totalX = float(skysize-1)
    domegrid = numpy.zeros([skysize, skysize], 'float')
    # Gradient is from top to bottom. 
    for i in range(0, len(pixgrid[0])):
        for j in range(0, len(pixgrid[1])):
            xfrac = j/totalX
            domegrid[i][j] = domeval + domeval*domegradient*xfrac
    pixgrid = pixgrid*domegrid #+ pixgrid*domeval
    return pixgrid, domegrid

def add_pixscale(pixgrid):
    # Add effects of pixel variance. 
    skygrid = numpy.zeros([skysize, skysize], 'float')
    pixcenX = skysize-1
    pixcenY = skysize-1
    totaldistance = numpy.sqrt(pixcenX**2+pixcenY**2)
    for i in range(0, len(pixgrid[0])):
        for j in range(0, len(pixgrid[1])):
            distance = numpy.sqrt((i-pixcenX)**2 + (j-pixcenY)**2)
            # Calculate pixel area, this is what skygrid will be - showing pixel area variation
            pixside = 1 + pixscale_change*(distance/totaldistance)
            skygrid[i][j] = pixside * pixside
    #normforsky = pixgrid.mean()
    pixgrid = pixgrid * skygrid
    #skygrid = skygrid * normforsky
    return pixgrid, skygrid
        
def add_ghost(pixgrid, rawdetector):
    # Add a ghost pupil image (for simplicity sake) in a ring, with value proportional to total counts.
    # should be added after dome screen variance, if applicable. 
    totalcounts = pixgrid.sum()
    ghostring = numpy.zeros([skysize, skysize], 'float')
    pixcenX = skysize-1
    pixcenY = skysize-1
    # Set pix limits for ghost
    ghostringdiam_outer = skysize/1.2
    ghostringdiam_inner = skysize/2.7
    ghostarea = (ghostringdiam_outer**2 - ghostringdiam_inner**2)*numpy.pi
    print "ghost: total %f, area %f, added %f per pix" %(totalcounts, ghostarea, totalcounts/ghostarea*ghostscale)
    # Generate ghost pupil image.
    for i in range(0, len(pixgrid[0])):
        for j in range(0, len(pixgrid[1])):
            distance = numpy.sqrt((i-pixcenX)**2 + (j-pixcenY)**2)
            if distance>ghostringdiam_outer:
                ghostring[i][j] = 0
            elif distance<ghostringdiam_inner:
                ghostring[i][j] = 0
            else:
                ghostring[i][j] = totalcounts/ghostarea*ghostscale*rawdetector[i][j]
    # Add ghost to pixel image.
    pixgrid = pixgrid + ghostring
    return pixgrid, ghostring
            
    

def make_flatfigure():
    domeval = 1
    detector = make_pixgrid()
    flat = numpy.copy(detector)
    flat, dome = add_domewvariance(flat, domeval)
    flat, pixscale = add_pixscale(flat)
    flat, ghost = add_ghost(flat, detector)
    print "created flat and added sky, dome, ghost"
    # Set up figure. 
    fig = pylab.figure()
    fig.subplots_adjust(wspace=0.45, hspace=0.1, top=0.93, bottom=0.1)
    ax = pylab.subplot(421)
    pylab.figtext(0.05, 0.83, "Detector")
    im = ax.imshow(detector, origin='lower')    
    nsteps = 2.0
    step = (detector.max()-detector.min())/nsteps
    ticks = numpy.arange(detector.min(), detector.max()+step/2.0, step)
    pylab.colorbar(im, shrink=0.75, aspect=10, ticks=ticks, format='%.3f')
    pylab.grid(which='major')
    ax.tick_params(axis='x', labelbottom='off')
    ax = pylab.subplot(423)
    pylab.figtext(0.05, 0.61, "Dome Var.")
    im = ax.imshow(dome, origin='lower')
    step = (dome.max()-dome.min())/nsteps
    ticks = numpy.arange(dome.min(), dome.max()+step/2.0, step)
    pylab.colorbar(im, shrink=0.75, aspect=10, ticks=ticks, format='%.3f')
    pylab.grid(which='major')
    ax.tick_params(axis='x', labelbottom='off')
    ax = pylab.subplot(425)
    pylab.figtext(0.05, 0.4, "Pixel Scale")
    im = ax.imshow(pixscale, origin='lower')
    step = (pixscale.max()-pixscale.min())/nsteps
    ticks = numpy.arange(pixscale.min(), pixscale.max()+step/2.0, step)
    pylab.colorbar(im, shrink=0.75, aspect=10, ticks=ticks, format='%.3f')
    pylab.grid(which='major')
    ax.tick_params(axis='x', labelbottom='off')
    ax =pylab.subplot(427)
    pylab.figtext(0.05, 0.2, "Ghosting")
    im = ax.imshow(ghost, origin='lower')
    step = (ghost.max()-ghost.min())/nsteps
    ticks = numpy.arange(ghost.min(), ghost.max()+step/2.0, step)
    pylab.colorbar(im, shrink=0.75, aspect=10, ticks=ticks, format='%.4f')
    pylab.grid(which='major')
    pylab.xticks(rotation=-45)
    ax = pylab.subplot(222)
    pylab.figtext(0.65, 0.92, "Dome Flat")
    im = ax.imshow(flat, origin='lower')
    nsteps = 4
    step = (flat.max()-flat.min())/nsteps
    ticks = numpy.arange(flat.min(), flat.max()+step/2.0, step)
    pylab.colorbar(im, shrink=0.8, ticks=ticks, format='%.3f')
    pylab.grid(which='major')
    ax = pylab.subplot(224)
    pylab.figtext(0.65, 0.48, "Illum. Corr.")
    ic = (flat - ghost)/dome/pixscale/flat
    im = ax.imshow(ic, origin='lower')
    step = (ic.max()-ic.min())/nsteps
    ticks = numpy.arange(ic.min(), ic.max()+step/2.0, step)
    pylab.colorbar(im, shrink=0.8, ticks=ticks, format='%.3f')
    pylab.grid(which='major')
    #return flat, ic

def make_skywflat():
    # make the flat field
    domeval = 10000
    detector = make_pixgrid()
    flat = numpy.copy(detector)
    flat, dome = add_domewvariance(flat, domeval)
    flat, pixscale = add_pixscale(flat)
    flat, ghost = add_ghost(flat, detector)
    ic = (flat - ghost)/dome/pixscale/flat
    ic = ic / ic.mean()
    flat = flat/flat.mean()
    # add stars to the 'sky'
    sky = numpy.copy(detector)
    star1 = [30, 30, 5000]
    star2 = [200, 200, 5000]
    sky = add_bkgd(sky, 1000)
    sky, pixscale2 = add_pixscale(sky)
    star = add_star(star1, detector)
    sky = sky + star
    star = add_star(star2, detector)
    sky = sky + star
    sky, ghost2 = add_ghost(sky, detector)
    # plot the results.
    fig = pylab.figure()
    fig.subplots_adjust(wspace=0.25, hspace=0.01, left=0.08, right=0.94)
    ax = pylab.subplot(231)
    pylab.title("Dome Flat, no I.C.", fontsize='medium')
    im = ax.imshow(flat, origin='lower')
    pylab.colorbar(im, shrink=0.65, format="%.3f")
    pylab.grid(which='major')
    pylab.xlim(0,skysize-1)
    pylab.ylim(0,skysize-1)
    ax.tick_params(axis='x', labelbottom='off')
    ax = pylab.subplot(232)
    pylab.title("Raw Sky", fontsize='medium')
    im = ax.imshow(sky, origin='lower')
    pixgrid = sky
    starval1, bkgd1 = measure_star(pixgrid, star1)
    starval2, bkgd2 = measure_star(pixgrid, star2)
    pylab.annotate("%d" %(round(starval1)), [star1[0]+5,star1[1]+12])
    pylab.annotate("%d" %(round(starval2)), [star2[0]-50,star2[1]+12])
    step = (sky.max() - sky.min())/5.0
    ticks = numpy.arange(sky.min()-step, sky.max()+2*step, step)
    ticks = numpy.round(ticks)
    pylab.colorbar(im, shrink=0.65, ticks=ticks, format='%d')
    pylab.grid(which='major')
    pylab.xlim(0,skysize-1)
    pylab.ylim(0,skysize-1)
    ax.tick_params(axis='x', labelbottom='off')
    ax.tick_params(axis='y', labelleft='off')
    ax = pylab.subplot(233)
    pylab.title("Sky/Dome, no IC", fontsize='medium')
    pixgrid = (sky-ghost2)/flat
    im = ax.imshow(pixgrid, origin='lower')
    step = (pixgrid.max() - pixgrid.min())/5.0
    ticks = numpy.arange(pixgrid.min()-step, pixgrid.max()+2*step, step)
    ticks = numpy.round(ticks)
    pylab.colorbar(im, shrink=0.65, ticks=ticks)
    starval1, bkgd1 = measure_star(pixgrid, star1)
    starval2, bkgd2 = measure_star(pixgrid, star2)
    pylab.annotate("%d" %(round(starval1)), [star1[0]+5,star1[1]+12])
    pylab.annotate("%d" %(round(starval2)), [star2[0]-50,star2[1]+12])
    ax.tick_params(axis='x', labelbottom='off')
    ax.tick_params(axis='y', labelleft='off')
    pylab.grid(which='major')
    pylab.xlim(0,skysize-1)
    pylab.ylim(0,skysize-1)
    ax = pylab.subplot(234)
    pylab.title("Dome Flat & IC", fontsize='medium')
    im = ax.imshow(flat*ic, origin='lower')
    pylab.colorbar(im, shrink=0.65, format='%.3f')
    pylab.grid(which='major')
    pylab.xlim(0,skysize-1)
    pylab.ylim(0,skysize-1)    
    pylab.xticks(rotation=-45)
    ax = pylab.subplot(235)
    pylab.title("Raw Sky", fontsize='medium')
    im = ax.imshow(sky, origin='lower')
    pixgrid = sky
    starval1, bkgd1 = measure_star(pixgrid, star1)
    starval2, bkgd2 = measure_star(pixgrid, star2)
    pylab.annotate("%d" %(round(starval1)), [star1[0]+5,star1[1]+12])
    pylab.annotate("%d" %(round(starval2)), [star2[0]-50,star2[1]+12])
    step = (sky.max() - sky.min())/5.0
    ticks = numpy.arange(sky.min()-step, sky.max()+2*step, step)
    ticks = numpy.round(ticks)
    pylab.colorbar(im, shrink=0.65, ticks=ticks)
    pylab.grid(which='major')
    pylab.xlim(0,skysize-1)
    pylab.ylim(0,skysize-1)
    ax.tick_params(axis='y', labelleft='off')
    pylab.xticks(rotation=-45)
    ax = pylab.subplot(236)
    pylab.title("Sky/Dome & IC", fontsize='medium')
    pixgrid = (sky-ghost2)/flat/ic        
    im = ax.imshow(pixgrid, origin='lower')
    step = (pixgrid.max() - pixgrid.min())/5.0
    ticks = numpy.arange(pixgrid.min()-step, pixgrid.max()+step*2, step)
    ticks = numpy.round(ticks)
    pylab.colorbar(im, shrink=0.65, ticks=ticks)
    pylab.grid(which='major')
    ax.tick_params(axis='y', labelleft='off')
    pylab.xticks(rotation=-45)
    starval1, bkgd1 = measure_star(pixgrid, star1)
    starval2, bkgd2 = measure_star(pixgrid, star2)
    pylab.annotate("%d" %(round(starval1)), [star1[0]+5,star1[1]+12])
    pylab.annotate("%d" %(round(starval2)), [star2[0]-50,star2[1]+12])
    pylab.xlim(0,skysize-1)
    pylab.ylim(0,skysize-1)

