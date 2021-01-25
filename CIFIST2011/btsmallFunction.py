import gzip
from numpy import loadtxt, arange, array, zeros, log, exp
from numpy import logical_or, isnan, amin, amax
from cPickle import load
from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import InterpolatedUnivariateSpline

infile = '/home/jgallimo/Documents/DREAMZ_RV_TESTS/CIFIST2011/btsmall.pkl.gz'
fin = gzip.open(infile)
protostarModels = load(fin)
fin.close()

protoW = protostarModels['wave']
tList = protostarModels['T']
gList = protostarModels['g']
protoModels = protostarModels['models']

tDex = InterpolatedUnivariateSpline(tList, arange(len(tList)), k=1)
gDex = InterpolatedUnivariateSpline(gList, arange(len(gList)), k=1)

indwave = arange(len(protoW))
iwf = InterpolatedUnivariateSpline(protoW, indwave, k=1)

def protoFunction(t, g):
    it = tDex(t)
    ig = gDex(g)
    temp = array([it, ig], dtype=float)
    coords = zeros((3, len(protoW)))
    coords[0,:] = indwave
    for j in arange(len(temp)):
        coords[j+1,:] = temp[j]
    result = map_coordinates(protoModels, coords, order=1, \
                             mode='constant', cval=0)
    return result
    
