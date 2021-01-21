#!/usr/bin/env python
import matplotlib as mpl
from scipy.interpolate import RegularGridInterpolator
mpl.use('Agg')
import matplotlib.pyplot as mp
mp.ioff()
mpl.rc('font',family='Times New Roman')
import time
import amoeba
from scipy.optimize import differential_evolution, fmin, basinhopping
from pyfits import getdata
from matplotlib.pyplot import clf, axes, savefig
from numpy import inf, loadtxt, exp, log, logical_and, sqrt, copy, sum
from numpy import median, array, zeros, ones, linspace, arange
from numpy import amin, amax, pi, ceil, isnan, polyfit, polyval
from numpy import isinf, diff, all, RankWarning, average, std, isfinite
from numpy import argmax, nan, append
from numpy.random import normal, uniform
from gzip import open as zopen
from dreamZPT.wrapper import DREAMZPT
#from dreamZPT.wrapper import DREAMZS
from astrolib.FindHelio import findCorr
from dreamZPT.diagnostics import autoburn
from cPickle import load, dump
from scipy.ndimage.interpolation import map_coordinates
#from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import InterpolatedUnivariateSpline
from glob import glob
from pylab import *
import gzip
import math
import os

import warnings
warnings.simplefilter('ignore', RankWarning)
from scipy.interpolate import BarycentricInterpolator

from scipy.signal import gaussian, fftconvolve
from scipy.ndimage.filters import gaussian_filter1d, convolve1d, \
    uniform_filter1d

modelList = sorted(glob('./KBand/*.gz'))
# debugging override
# modelList = ['/home/jgallimo/Documents/dreamzBrownDwarfs/KBand/lte014.5-4.5-kband.txt.gz']

# load a token xmodel, ymodel
xmodel, ymodel = loadtxt(modelList[0], unpack=True)

import argparse
parser = argparse.ArgumentParser(description='Fit a Brown Dwarf Spectrum')
parser.add_argument('fileName', help='File Name')
args = parser.parse_args()
fileName = args.fileName

#fileName = './Data/PSOJ318_comb.fits'
niter = 100000
sysUnc = 0.
#fixLSF = False
thin = 10
ncr = 3
npoly = 2
navg = 4
noseed = True
nchains = 3
modelfileName = './CIFIST2011/btsettl-3.pkl.gz'


def reportElapsedTime(startTime):
    elapsed = time.time() - startTime
    rElapsed = elapsed
    ts = 'seconds'
    if elapsed > 60.:
        rElapsed = elapsed / 60.
        ts = 'minutes'
    if elapsed > 3600:
        rElapsed = elapsed / 3600.
        ts = 'hours'
    if elapsed > 3600*24.:
        rElapsed = elapsed / 3600. / 24.
        ts = 'days'
    if elapsed > 3600.*24.*365.:
        rElapsed = elapsed / (3600.*24*365.)
        ts = 'years'
    print 'Elapsed time: ', round_sigfigs(rElapsed, 3), ts
    return


fin = gzip.open(modelfileName)
protostarModels = load(fin)
fin.close()

protoW = protostarModels['wave']

#fileName2 = fileName1 
#fileName = 'PSOJ318_comb'

seed = None # './PSOJ318_seed.txt'

bestscore = -inf
bestmodel = None

thismod = modelfileName.split('/')
thismod = thismod[-1].split('.')
suffix = 'dreamZS-' + thismod[0]

debug = False
c = 299792.458 # km/s
MAXVEL = 300. # km/s
OVRSAMPL = 10.
if '2015' in modelfileName:
    OVRSAMPL=5.
#bandLow = 2.1
#bandHigh = 2.40

'''
Load Data
'''
aStarData, hdr = getdata(fileName, header=True)
RP = hdr['RP']
order = hdr['ORDERS'].replace(" ","")
if "_" in fileName:
    sourceName = hdr['OBJECT'].replace(" ","")
#sourceName =fileName.split("_")[0]
else:
    sourceName = fileName

outdir = sourceName + '-Results'
if not os.path.exists(outdir):
    os.makedirs(outdir)


xStar = aStarData[0,:]
dx = median(xStar[1] - xStar[0])
yStar = aStarData[1,:]
errStar = aStarData[2,:]

# override for PSOJ
'''
idx = (xStar >= 2.27) * (xStar <= 2.33)
xStar = xStar[idx]
yStar = yStar[idx]
errStar = errStar[idx]
'''

nGrid = len(yStar)
nominalFWHM = median(xStar) / RP


sigh = isnan(yStar)
yStar[sigh] = -9999.

good = isfinite(yStar) * (yStar > 0) * isfinite(errStar)
bad = ~good

#initialize wavelengths for interpolation
wmin = amin(xStar) - 0.01 #minimum wavelength
wmax = amax(xStar) + 0.01 #maximum wavelength
if wmax > amax(protoW):
    wmax = amax(protoW)
if wmin < amin(protoW):
    wmin = amin(protoW)

modelWavelength = logical_and(protoW>= wmin, protoW<=wmax)
protoW = protoW[modelWavelength]
mfwSub = linspace(amin(protoW), amax(protoW), OVRSAMPL*len(protoW))
dxos = mfwSub[1] - mfwSub[0]

tList = protostarModels['T']
gList = protostarModels['g']

protoModels = protostarModels['models']
protoModels = protoModels[modelWavelength, :, :] # log models

tDex = InterpolatedUnivariateSpline(tList, arange(len(tList)), k=1)
gDex = InterpolatedUnivariateSpline(gList, arange(len(gList)), k=1)

indwave = arange(len(protoW))
iwf = InterpolatedUnivariateSpline(protoW, indwave, k=1)

# load transmission curve
xtransg, ytransg = loadtxt('transmissionCurve.txt', unpack=True)
dltransg = xtransg[1] - xtransg[0]
yyg = InterpolatedUnivariateSpline(xtransg, ytransg, k=3)(mfwSub)
yyg[yyg < 0] = 0.0

# convenience: calculate velocity grid outside of model function
l0 = log(amin(protoW))
l1 = log(amax(protoW))
# 10 times oversampled
mwlVgrid = exp(linspace(l0, l1, len(protoW)*OVRSAMPL))
vpix = c * median((mwlVgrid[1:] - mwlVgrid[:-1]) / mwlVgrid[1:])


# note: I'm finding residual scatter of ~ 10%
# I'll add it in quadrature
eStar = sqrt((sysUnc*yStar)**2 + errStar**2)
vhelio = findCorr(fileName)

def rolling_mean(x, n):
    return fftconvolve(x, ones((n,))/float(n), mode='same')

def robustpolyfit(x, y, deg, rcond = None, full=False, w=None, cov=False, niter=3, nrej=4):
    '''
    niter = number of rejection iterations
    nrej = number of sigma rejection
    '''
    '''
    perform initial polyfit
    '''
    n = len(x)
    if w == None:
        weight = ones(n)
    else:
        weight = w
    tgood = ones(n, dtype='bool')
    tgood[~isfinite(y) + ~isfinite(weight)] = False
    tgood[y < 0] = False
    c = polyfit(x[tgood], y[tgood], deg, rcond = rcond, w=weight[tgood], cov=cov)
    m = polyval(c, x)
    r = weight*(y - m)
    sig = std(r[tgood])
    iter = 0
    tgood[abs(r) > nrej*sig] = False
    newReject = n - sum(tgood)
    while(iter < niter and newReject > 0):
        c = polyfit(x[tgood], y[tgood], deg, \
                    rcond = rcond, w=weight[tgood], cov=cov)
        m = polyval(c, x)
        r = weight*(y - m)
        sig = std(r[tgood])
        tgood[abs(r) > nrej*sig] = False
        reject = n - sum(tgood)
        newReject = reject - newReject
        iter += 1
    if w != None:
        c = polyfit(x[tgood], y[tgood], deg, rcond = rcond, w=weight[tgood], cov=cov, full=full)
    else:
        c = polyfit(x[tgood], y[tgood], deg, rcond = rcond, w=weight[tgood], cov=cov, full=full)    
    return c

def round_sigfigs(num, sig_figs):
    """Round to specified number of sigfigs.

    >>> round_sigfigs(0, sig_figs=4)
    0
    >>> int(round_sigfigs(12345, sig_figs=2))
    12000
    >>> int(round_sigfigs(-12345, sig_figs=2))
    -12000
    >>> int(round_sigfigs(1, sig_figs=2))
    1
    >>> '{0:.3}'.format(round_sigfigs(3.1415, sig_figs=2))
    '3.1'
    >>> '{0:.3}'.format(round_sigfigs(-3.1415, sig_figs=2))
    '-3.1'
    >>> '{0:.5}'.format(round_sigfigs(0.00098765, sig_figs=2))
    '0.00099'
    >>> '{0:.6}'.format(round_sigfigs(0.00098765, sig_figs=3))
    '0.000988'
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0



def lsf_rotate(deltav, vsini, epsilon=0.6):  
    e1 = 2.0 * (1.0 - epsilon)
    e2 = pi * epsilon / 2.0
    e3 = pi * (1.0 - epsilon / 3.0)
    npts = int(ceil(2.0 * vsini / deltav))
    if npts % 2 == 0: 
        npts += 1
    nwid = npts / 2
    xx = arange(npts) - 1.0*nwid
    xx *= deltav / vsini
    velgrid = xx * vsini
    x1 = abs(1.0 - xx*xx)
    lsf = (e1*sqrt(x1) + e2 * x1)/e3
    return velgrid, lsf / sum(lsf) # experiment in normalization

#protoFunction = RegularGridInterpolator((protoW,tList, gList), protoModels, method='linear')


def protoFunction(t, g):
    it = tDex(t)
    ig = gDex(g)
    temp = array([it, ig], dtype=float)
    coords = zeros((3, len(protoW)))
    coords[0,:] = indwave
    for j in arange(len(temp)):
        coords[j+1,:] = temp[j]
    result = map_coordinates(protoModels, coords, order=1, mode='constant', cval=0)
    return result # log models


def strictly_increasing(xx, thresh):
    dx = diff(xx)
    return all(dx > thresh)

def weighAnchor(nGrid, anchors):
    xx = array(linspace(0.,1.,npoly+1)*(nGrid-1))#.astype(int)
    return xx

def xgridFromAnchors(nGrid, p):
    anchors = p[6:]
    pivots = weighAnchor(nGrid, anchors)
    x = BarycentricInterpolator(pivots,anchors)(arange(nGrid))
    return x


def logPrior(pp): # p is an array of parameters
    p = zeros(7+npoly)
    p[0:4] = pp[0:4]
    p[4] = t
    p[5] = g
    p[6:] = pp[4:]
    fwhm = p[0]*1000.
    alpha = p[1]
    vsini = p[2]
    rvmag = abs(p[3])
    penalty = 0.0
    if not strictly_increasing(p[6:], 0.000):
        return -inf
    if fwhm < 0:
        return -inf
    if fwhm > 0.3:
        penalty -= (fwhm - 0.3)**2 / (0.01)**2 # roughly the expected fwhm
    if alpha > 1.0 or alpha < 0.: # that shouldn't happen! Or should it?
        return -inf
    if vsini < 0.:
        return -inf
    if vsini > 100.:
        penalty -= (vsini - 100)**2 / (1.)**2
    if rvmag > 100.: # avoid absurdly high rvs
        penalty -= (rvmag - 100.)**2 / (1.)**2
    if p[6] <= 1.8: # keep things in the K-band
        return -inf
    if p[-1] >= 2.5:
        return -inf
    
    return penalty 


def negchisq(x, data=(xmodel, ymodel)):
    logp = logPrior(x)
    if isinf(logp):
        return -inf
    logl = modelLogP(x, data[0], data[1], y=yStar, e=eStar, sysUnc=sysUnc) + logp
    chisq = 2. * logl
    return chisq 

def chisq(x, xmodel, ymodel):
    #logp = logPrior(x)
    #if isinf(logp):
    #    return inf
    logl = modelLogP(x, xmodel, ymodel, y=yStar, e=eStar, sysUnc=sysUnc)# + logp
    chisq = -2. * logl
    #print chisq, x
    return chisq 

def model(pp, xmodel, ymodel, y=yStar, getPolyCoef=False, full=False):
    # calculate wavelength grid from parameters
    # (refine wavelength calibration)
    p = zeros(7+npoly)
    p[0:4] = pp[0:4]
    p[4] = t
    p[5] = g
    p[6:] = pp[4:]
    x = xgridFromAnchors(nGrid, p)
    if not strictly_increasing(x, 0.):
        model2 = 0.0 * x
        if getPolyCoef:
            coef = zeros(4)
            return model2, coef
        return model2

    fwhm = p[0]
    alpha = p[1]
    vsini = p[2]
    rv = p[3]

    #mfl = protoFunction(p[4], p[5])
    mfl = ymodel
    protoW = xmodel
    mflVgrid = InterpolatedUnivariateSpline(protoW, mfl, k=3)(mwlVgrid)
    '''
    mflVgrid = mwlVgrid * 0.0
    for i in range(len(mwlVgrid)):
        mflVgrid[i] = protoFunction([mwlVgrid[i],p[4], p[5]])
    '''
    if vsini > vpix:
        velgrid, vsiniSmooth = lsf_rotate(vpix, vsini)
        mflVgridC = fftconvolve(mflVgrid, vsiniSmooth, mode='same')
        #mflVgridC = convolve1d(mflVgrid, vsiniSmooth, mode='nearest')
    else:
        mflVgridC = copy(mflVgrid)

    # shift back to model grid, observer's reference frame
    mwlVgridS = mwlVgrid * (1.0 + (rv-vhelio) / c) 
    mflS = InterpolatedUnivariateSpline(mwlVgridS, mflVgridC)(mfwSub)

    # Resample sky transmission to model wavelength grid
    stp1 = yyg**alpha
    #print 'Max stp1 = ', amax(stp1), alpha
    # apply sky transmission
    model0 = stp1 * mflS

    # Smooth by lsf
    sigma = fwhm / dxos / 2.355 # convert FWHM to sigma
    kernel = gaussian(int(10*sigma+1), sigma)
    kernel /= sum(kernel)
    model1 = fftconvolve(model0, kernel, mode='same') 
    #model1 = gaussian_filter1d(model0, sigma)

    # interpolate back to data grid
    try: # check if the wavelength solution is running away
        model2 = InterpolatedUnivariateSpline(mfwSub, model1, k=2)(x) 
    except:
        model2 = 0.0 * x
        if getPolyCoef:
            coef = zeros(4)
            return model2, coef
        return model2
    
    # determine normalization, scale protostar model to data
    # Apply a 3rd order polynomial flux scaling
    j = arange(len(yStar)) / float(len(yStar) - 1)

    # Amoeba technique
    '''
    try:

    coef = robustpolyfit(j[good], yStar[good]/model2[good], 3, \
    w=errStar[good]/model2[good])
    model3 = polyval(coef, j)*model2

    except:
    '''
    coef = polyfit(j[good], rolling_mean(yStar[good]/model2[good], 20), 3)
    model3 = polyval(coef, j)*model2
    scl = median(yStar[good])/median(model3[good])
    model3 *= scl

    if getPolyCoef:
        return model3, coef
    if full:
        # calculate smoothed transmission spectrum
        smootrans = fftconvolve(stp1, kernel, mode='same') 
        #smootrans = convolve1d(stp1, kernel, mode='nearest')
        smootrans = InterpolatedUnivariateSpline(mfwSub, smootrans, k=2)(x)
        # calculate smoothed star-only spectrum, without telluric features
        model0 = fftconvolve(mflS, kernel, mode='same')
        #model0 = convolve1d(mflS, kernel, mode='nearest')
        model2 = polyval(coef, j)*InterpolatedUnivariateSpline(mfwSub, \
                                                 mflS, k=2)(x)         
        '''
        model2 = polyval(coef, model2)*InterpolatedUnivariateSpline(mfwSub, \
                                                 mflS, k=2)(x)         
        '''
        return model3, coef, smootrans, model2
    return model3


def modelLogP(p, xmodel, ymodel, y=yStar, e=eStar, sysUnc=sysUnc, good=good, data=None):
    global bestscore, bestmodel
    m = model(p, xmodel, ymodel)
    z = (y-m) / e
    #fracErr =(y-m)/m
    #idx = abs(fracErr) < sysUnc 
    #z[idx] = sqrt(1.0 - 9. / len(yStar))
    logp = -sum(z[good]*z[good])/2.0
    if isinf(logp) or isnan(logp): logp = -inf
    if logp > bestscore:
        k = len(fitvector)
        bestscore = logp
        bestmodel = copy(m)

        print '---------------'
        print 'New optimum'
        print 'chisq: ', round_sigfigs(-2.0*bestscore, 4)
        print '  dof: ', round_sigfigs(len(yStar) - k, 4)
        print ' fwhm = ', round_sigfigs(p[0]*1000., 4), 'nm'
        print '   RV = ', round_sigfigs(p[3], 4), 'km/s'
        print 'vsini = ', round_sigfigs(p[2], 4), 'km/s'
        print 'alpha = ', round_sigfigs(p[1], 4)
        print 'lambda = ', p[4:]

    #print logp, p
    #print -2.0 * logp, p
    return logp


#Initialize Z array
# Require a seed for this version.
# Seed from dreamBDopt.py
nn =  7 + npoly # number of parameters
fitvector = arange(nn)

outname = outdir + '/' + sourceName + '-' + order
fp = open(outname + '-gtmap2.txt', 'w')

for m in modelList:
    xmodel, ymodel = loadtxt(m, unpack=True)
    t = float(m.split('/')[-1][3:8])
    g = float(m.split('/')[-1][9:12])

    #p0 = [0.00010, 0.869, 5.0, 35.0, xStar[0], median(xStar), xStar[-1]]
    # Here's where I put in initial guesses
    # lsf, tau, rv, vsini, (wavelengths for parabolic corrections)
    p0 = [0.00010, 0.5, 0.0, 25.0, xStar[0], median(xStar), xStar[-1]]
    scl=[0.2*p0[0], 0.1*p0[1], 15.0, 25.0, p0[0]*2.0, p0[0]*2.0, p0[0]*2.0]

    aresult = amoeba.amoeba(p0, scl, negchisq, data=(xmodel, ymodel))
    p0 = aresult[0]
    scl = list(array(scl)/2.0)
    aresult = amoeba.amoeba(p0, scl, negchisq, data=(xmodel, ymodel))
    p0 = aresult[0]

    print -aresult[1],
    print t, g,
    for thing in aresult[0]:
        print thing,
    print ''
    '''
    print 'Starting Basin-Hopping'
    kwargs = {}
    kwargs['args'] = (xmodel, ymodel)
    kwargs['method'] = 'Nelder-Mead'
    bresult = basinhopping(chisq, p0, minimizer_kwargs=kwargs)
    print  bresult['fun'],
    print  t, g,
    for thing in bresult['x']:
        print  thing,
    print  ''
    '''
    if t >= 12. and t <= 28. and g >= 3.5 and g <= 5.5:
        print >>fp, -aresult[1],
        print >>fp, t, g,
        for thing in aresult[0]:
            print >>fp, thing,
        print >>fp, ''
        '''
        print  >>fp, bresult['fun'],
        print  >>fp, t, g,
        for thing in bresult['x']:
            print  >>fp, thing,
        print  >>fp, ''
        '''

    '''
    print bresult['fun'],
    print t, g,
    for thing in bresult['x']:
        print thing,
    print ''
    '''

fp.close()
'''
# Code to assess individual models
for t in [14.0, 14.5, 15.0]:
    g = 4.5
    m = './KBand/lte%05.1f-%3.1f-kband.txt.gz' % (t, g)
    xmodel, ymodel = loadtxt(m, unpack=True)
    p0 = [0.161e-3, 0.841, 15.0, -6.0, 2.26783, 2.299085, 2.3305]
    scl=[0.2*p0[0], 0.1*p0[1], 15.0, 25.0, p0[0]*2.0, p0[0]*2.0, p0[0]*2.0]
    aresult = amoeba.amoeba(p0, scl, negchisq, data=(xmodel, ymodel))
    p0 = aresult[0]
    scl = list(array(scl)/2.0)
    aresult = amoeba.amoeba(p0, scl, negchisq, data=(xmodel, ymodel))
    p0 = aresult[0]
    pp = zeros(7 + npoly)
    pp[0:4] = p0[0:4]
    pp[4] = t
    pp[5] = g
    pp[6:] = p0[4:]
    xx = xgridFromAnchors(nGrid, pp)
    modl = model(p0, xmodel, ymodel)
    clf()
    plot(xx, yStar, 'k-')
    plot(xx, modl, 'r-')
    xlabel('Wavelength (microns)')
    ylabel('Normalized Flux Density')
    savefig(outname + '-amoeba-model-%05.1f-%3.1f.pdf' % (t,g))
    fp = open(outname + '-amoeba-model-%05.1f-%3.1f.txt' % (t,g), 'w')
    print >>fp, '# lambda(corrected) model data unc'
    print >>fp, '# fwhm = ', pp[0]
    print >>fp, '# tau = ', pp[1]
    print >>fp, '# vsini = ', pp[2]
    print >>fp, '# vr = ', pp[3]
    print >>fp, '# wavelengths = ', pp[6:]
    print >>fp, '# chisq = ', -aresult[1]

    for i in range(len(xx)):
        print >>fp, xx[i], modl[i], yStar[i], errStar[i]
    fp.close()
'''
