#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as mp
mp.ioff()
mpl.rc('font',family='Times New Roman')
#from pyfits import getdata
from astropy.io.fits import getdata
from matplotlib.pyplot import clf, axes, savefig
from numpy import inf, loadtxt, exp, log, logical_and, sqrt, copy, sum
from numpy import median, array, zeros, ones, linspace, arange
from numpy import amin, amax, pi, ceil, isnan, polyfit, polyval
from numpy import isinf, diff, all, RankWarning, average, std, isfinite
from numpy import argmax, nan, append
from numpy.random import normal, uniform
from gzip import open as zopen
from dreamZPT.wrapper import DREAMZPT
from dreamZPT.wrapper import DREAMZS
from astrolib.FindHelio import findCorr
from dreamZPT.diagnostics import autoburn
#from cPickle import load, dump
import pickle
from scipy.ndimage.interpolation import map_coordinates
#from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import InterpolatedUnivariateSpline

# Example run:
# ./dreamzsBD8.py ./Data/VB10_order32.fits --niter 30000  --nchains 3


import gzip
import math
import os

import warnings
warnings.simplefilter('ignore', RankWarning)
from scipy.interpolate import BarycentricInterpolator

from scipy.signal import gaussian, fftconvolve


import argparse
parser = argparse.ArgumentParser(description='Fit a Brown Dwarf Spectrum')
parser.add_argument('fileName', help='File Name')
#parser.add_argument('--fixLSF', help='Don\'t fit the LSF', default=False, \
#                    action='store_true')
parser.add_argument('--modelfile', help='Name of model grid (default: ./CIFIST2011/btsettl-4.pkl.gz', default='./CIFIST2011/btsettl-4.pkl.gz')
parser.add_argument('--pton', help='Turn off parallel tempering? Default=False', default=False, action='store_true')
parser.add_argument('--niter', help='Number of iterations (default: 100000)', \
                   default=100000, type=int)
parser.add_argument('--nchains', \
                    help='DREAMZ # of chains(default: 3)', 
                    default=3, type=int)
parser.add_argument('--fsys', \
                   help='Systematic Fractional Uncertainty (default: 0)', 
                   default=0, type=float)
parser.add_argument('--thin', \
                    help='MCMC thinning interval (default: 5)', 
                    default=5, type=int)
parser.add_argument('--ncr', \
                    help='DREAMZ crossover rate (default: 3; full crossover = -1)', 
                    default=3, type=int)
parser.add_argument('--navg', \
                    help='DREAMZ # of chains to average for DE (default: 4)', 
                    default=4, type=int)
parser.add_argument('--tlow', help='Lower limit on surface tempeature (default = 850 K)', default=850., type=float)
parser.add_argument('--thigh', help='Lower limit on surface tempeature (default = 2800 K)', default=2800., type=float)


args = parser.parse_args()
fileName = args.fileName
niter = args.niter
sysUnc = args.fsys
tlow = args.tlow / 100.
thigh = args.thigh / 100.
#fixLSF = args.fixLSF
thin = args.thin
ncr = args.ncr
ptoff = args.pton == False
navg = args.navg
nchains = args.nchains
modelfileName = args.modelfile
'''
fileName = './Data/vb10_order37.fits'
niter = 10000
ptoff = True
sysUnc = 0.1
#fixLSF = False
thin = 10
ncr = 3
navg = 4
nchains = 5
modelfileName = './CIFIST2011/btsettl-3.pkl.gz'
'''
fin = gzip.open(modelfileName)
protostarModels = pickle.load(fin)
fin.close()

protoW = protostarModels['wave']

#fileName2 = fileName1 
#fileName = 'PSOJ318_comb'

bestscore = -inf
bestmodel = None

thismod = modelfileName.split('/')
thismod = thismod[-1].split('.')
sfx = ''
if sysUnc > 0: sfx = 'fsys-'
if ptoff:
    suffix = 'dreamZS-' + sfx + thismod[0]
else:
    suffix = 'dreamZPT-' + sfx + thismod[0]
print suffix

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
if "NIRSPEC" in hdr['INSTR']:
    sourceName = hdr['TARGNAME'].replace(" ","")
#sourceName =fileName.split("_")[0]
else:
    sourceName = hdr['OBJECT'].replace(" ","")

xStar = aStarData[0,:]
dx = median(xStar[1] - xStar[0])
yStar = aStarData[1,:]
errStar = aStarData[2,:]
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

# temperature filter
idx = (tList >= tlow) * (tList <= thigh)
tList = tList[idx]
protoModels = protoModels[:, idx, :]
trange = (min(tList), max(tList))
grange = (min(gList), max(gList))

tDex = InterpolatedUnivariateSpline(tList, arange(len(tList)), k=1)
gDex = InterpolatedUnivariateSpline(gList, arange(len(gList)), k=1)

indwave = arange(len(protoW))
iwf = InterpolatedUnivariateSpline(protoW, indwave, k=1)

# load transmission curve
xtransg, ytransg = loadtxt('transmissionCurve_atrans.txt', unpack=True)
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
print 'HELIO VELOCITY:', vhelio

#vhelio = 28.1

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
    x=array(x)
    y=array(y)
    n = len(x)
    if w == None:
        weight = ones(n)
    else:
        weight = w
    tgood = ones(n, dtype='bool')
    tgood[~isfinite(y) + ~isfinite(weight)] = False
    #tgood[y < 0] = bool("False")
    #tgood[y < 0] = False
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

#protoFunction = RegularGridInterpolator((protoW,tList, gList), protoModels)


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


def logPrior(p): # p is an array of parameters
    fwhm = p[0]*1000.
    alpha = p[1]
    vsini = p[2]
    rvmag = abs(p[3])
    penalty = 0.0
    if not strictly_increasing(p[6:], 0.00):
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
    if p[4] < trange[0] or  p[4] > trange[1]: # temperatureGrid
        return -inf
    if p[5] < grange[0] or p[5] > grange[1]: # surfacegravityGrid
        return -inf
    if p[6] <= 1.8: # keep things in the K-band
        return -inf
    if p[-1] >= 2.5:
        return -inf

    return penalty 

def chisq(p):
    logp = logPrior(p)
    if isinf(logp):
        return inf
    logl = modelLogP(p, y=yStar, e=eStar, sysUnc=sysUnc) + logp
    chisq = -2. * logl
    return chisq 

def model(p, y=yStar, getPolyCoef=False, full=False):
    # calculate wavelength grid from parameters
    # (refine wavelength calibration)

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

    # get the model stellar atmosphere spectrum
    mfl = protoFunction(p[4], p[5])
    mflVgrid = InterpolatedUnivariateSpline(protoW, mfl, k=3)(mwlVgrid)
    if vsini > vpix:
        velgrid, vsiniSmooth = lsf_rotate(vpix, vsini)
        mflVgridC = fftconvolve(mflVgrid, vsiniSmooth, mode='same')
    else:
        mflVgridC = copy(mflVgrid)

    # shift back to model grid, observer's reference frame
    mwlVgridS = mwlVgrid * (1.0 + (rv-vhelio) / c) 
    mflS = InterpolatedUnivariateSpline(mwlVgridS, mflVgridC, k=3)(mfwSub)

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
    # Robust polyfit can be unstable. Here's the punt.
    try:
        coef = robustpolyfit(j, yStar/model2, 3, w=errStar/model2)
        model3 = polyval(coef, j)*model2
    except:
        yy = yStar / model2
        ok = isfinite(yy)
        coef = polyfit(j[ok], rolling_mean(yy[ok], 20), 3)
        model3 = polyval(coef, j)*model2
        scl = median(yStar[ok])/median(model3[ok])
        model3 *= scl
    if getPolyCoef:
        return model3, coef
    if full:
        # calculate smoothed transmission spectrum
        smootrans = fftconvolve(stp1, kernel, mode='same') 
        smootrans = InterpolatedUnivariateSpline(mfwSub, smootrans, k=2)(x)
        # calculate smoothed star-only spectrum, without telluric features
        model0 = fftconvolve(mflS, kernel, mode='same')
        model2 = polyval(coef, j)*InterpolatedUnivariateSpline(mfwSub, \
                                                 mflS, k=2)(x)         
        '''
        model2 = polyval(coef, model2)*InterpolatedUnivariateSpline(mfwSub, \
                                                 mflS, k=2)(x)         
        '''
        return model3, coef, smootrans, model2
    return model3


def modelLogP(p, y=yStar, e=eStar, sysUnc=sysUnc, good=good):
    global bestscore, bestmodel
    m = model(p)
    z = (y-m) / e
    '''
    if sysUnc > 0.:
        fracErr =(y-m)/m
        idx = abs(fracErr) < sysUnc 
        z[idx] = sqrt(1.0 - (1.0*len(fitvector)) / (1.0*len(yStar)))
    '''
    logp = -sum(z[good]*z[good])/2.0
    if isinf(logp) or isnan(logp): logp = -inf
    if logp > bestscore:
        k = len(fitvector)
        bestscore = logp
        bestmodel = m

        print '---------------'
        print 'New optimum'
        print 'chisq: ', round_sigfigs(-2.0*bestscore, 4)
        print '  dof: ', round_sigfigs(len(yStar) - k, 4)
        print ' fwhm = ', round_sigfigs(p[0]*1000., 4), 'nm'
        print '   RV = ', round_sigfigs(p[3], 4), 'km/s'
        print 'vsini = ', round_sigfigs(p[2], 4), 'km/s'
        print '    T = ', round_sigfigs(p[4]*100., 4), 'K'
        print ' logg = ', round_sigfigs(p[5], 4)
        print 'alpha = ', round_sigfigs(p[1], 4)
        print 'wavelength = ', p[6:]

    #print logp, p
    #print -2.0 * logp, p
    return logp


#Initialize Z array

outdir = sourceName + '-Results'

eps = 1.e-6

# Try initializing from gtmap
infile = sourceName + '-' + order + '-gtmap3.txt'
fp = open(outdir + '/' + infile, 'r')
init = fp.readlines()
fp.close()
ninit = len(init)
temp = init[0].split()[1:]
ndim = len(temp)
Z = zeros((ninit, ndim))
fitvector = arange(ndim)
npoly = ndim - 7

for i in range(ninit):
    temp = init[i].split()[1:]
    # messy because the order is wrong!
    Z[i,0] = float(temp[2])
    Z[i,1] = float(temp[3])
    Z[i,2] = float(temp[4])
    Z[i,3] = float(temp[5])
    Z[i,4] = float(temp[0])
    Z[i,5] = float(temp[1])
    for j in range(6, len(temp)):
        Z[i,j] = float(temp[j])
    print Z[i]
    

'''
for i in range(ninit):
    Z[i] = plows + uniform()*(phighs - plows)
    Z[i,0] = nominalFWHM * (1.0 + normal()*0.1)
    z = -1.0
    while z < 0:
        z = 0.7 + normal() * 0.1
    Z[i,1] = z
    Z[i,2] = abs(normal()*10.)
    Z[i,3] = normal()*10.
    print Z[i]
fitvector = arange(len(Z[0]))[~(plows == phighs)] 
'''
print fitvector
if ptoff:
    Z, ZP, successVec = DREAMZS(Z, modelLogP, logPrior, nchains, niter, \
                                thin=thin, ncr=ncr, navg = navg, \
                                args=(yStar, eStar, sysUnc), fitVector=fitvector)
else:
    XT = sqrt(2.0)**arange(nchains) # debug
    Z, ZP, ZT, successVec = DREAMZPT(Z, modelLogP, logPrior, XT, niter, \
                                     thin=thin, ncr=ncr, fitVector=fitvector, \
                                     navg = navg, args=(yStar, eStar, sysUnc), )

    # lose the high temperature solutions
    idx = ZT == 1
    Z = Z[idx]
    ZP = ZP[idx]

'''
Big change: save results to a Pickle file.
'''

idx = argmax(ZP)
p = Z[idx]
mdl, coef, trans, bdspec = model(p, full=True)
x = xgridFromAnchors(len(xStar), p)

db = {}
db['modelgrid'] = modelfileName
db['npoly'] = npoly
db['FWHM'] = abs(Z[:,0])
db['alpha'] = abs(Z[:,1])
db['vsini'] = abs(Z[:,2])
db['vr'] = Z[:,3]
db['T'] = Z[:,4] * 100.
db['logg'] = Z[:,5]
db['l0'] = Z[:,6]
db['l1'] = Z[:,7]
db['l2'] = Z[:,8]
chisq = -2.0*ZP
db['chisq'] = chisq
dof = len(yStar) - len(fitvector)
db['redChisq'] = amin(chisq) / float(dof)
db['dof'] = dof
db['burn'] = autoburn(chisq)
# re-flag data that we had interpolated over
yStar[bad] = nan
errStar[bad] = nan
db['errStar'] = errStar
db['xStar'] = xStar
db['wave'] = x
db['yStar'] = yStar
db['model'] = mdl
db['transmission'] = trans
db['BDmodel'] = bdspec
db['sysUnc'] = sysUnc
db['Z'] = Z
db['ZP'] = ZP

fp = zopen(outdir + '/' + sourceName + '-' + order + '-' + suffix + \
           '-Results.pkl.gz', 'w')
pickle.dump(db, fp)
fp.close()

# plot the model
p = Z[-1]
print 'p = ', p
x = xgridFromAnchors(len(yStar), p)
m = model(p)
clf()
left, width = 0.1, 0.8
bottom0, height0 = 0.1, 0.3
bottom1, height1 = bottom0+height0, 0.9 - bottom0-height0

spect = axes([left, bottom1, width, height1])
resid = axes([left, bottom0, width, height0])

resid.plot(x, yStar - m,'.')
resid.set_xlabel('Wavelength ($\\mu$m)')
resid.set_ylabel('Residuals')
#resid.set_yticks([-40,-20,0,20,40])
g = resid.get_xticklabels()

spect.plot(x, yStar, 'k')
spect.plot(x, m, 'r')
spect.set_xticklabels(g,visible=False)
spect.set_ylabel('Normalized $F_\\lambda$')
#spect.set_yticks(arange(400,800,50))

savefig(outdir + '/' + sourceName + '-' + order + '-' + suffix + \
        '-ModelFit.pdf')

# generate a wavelength-calibrated spectrum
# compile wavelength grids
#from dreamZPT.diagnostics import autoburn
#burn = autoburn(ZP)
#ZP = ZP[burn:]
#Z = Z[burn:]
nsamp = len(Z)
xgrid = zeros((nsamp, len(x)))
for i in range(nsamp):
    p = Z[i]
    x = xgridFromAnchors(len(yStar), p)
    xgrid[i] = x

meanXgrid = average(xgrid, axis=0)
stdXgrid = std(xgrid, axis=0)

fp = open(outdir + '/' + sourceName + '-' + order + '-' + suffix + \
          '_calibratedSpectrum.txt', 'w')
print >>fp, '# wavelength, error, flux, error'
for i in range(len(x)):
    print >>fp, meanXgrid[i], stdXgrid[i], yStar[i], errStar[i]
fp.close()

