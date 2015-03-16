import numpy as np
from scipy.interpolate import interp1d

def gen(wlMin, wlMax, dwl):
    
    errorLevel = np.array([[300.0, 0.01], [400.0, 0.002], [470.0, 0.001], [950.0, 0.001], [1200, 0.02]])

    (nwl, nwid) = errorLevel.shape

    errorSamp = np.zeros((nwl, nwid))
    errorSamp[:,0] = errorLevel[:,0]
    errorSamp[:,1] = np.random.randn(nwl) * errorLevel[:,1]

    wlOut = np.arange(wlMin, wlMax, dwl)
    errorOut = np.zeros(len(wlOut))

    interpFunc = interp1d(errorSamp[:,0], errorSamp[:,1], kind='cubic')

    errorOut = interpFunc(wlOut)

    responseArray = np.zeros((len(wlOut),2))
    responseArray[:,0] = wlOut
    responseArray[:,1] = errorOut + 1

    return errorSamp, responseArray

    
