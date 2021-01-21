#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as mp
mp.ioff()
from pylab import *
mpl.rc('font',family='Times New Roman')
# run in pylab environment
from scipy.stats import scoreatpercentile
from astrolib.FindHelio import findCorr
from dreamZPT.diagnostics import autoburn
import argparse
from cPickle import load
from gzip import open as zopen
from pyfits import getdata
import matplotlib
import os
import sys

parser = argparse.ArgumentParser(description='Print out the burned MCMC results in a text file')
parser.add_argument('pklfileName', help='File Name of the .pkl.gz file containing MCMC output')

args = parser.parse_args()
pklfileName = args.pklfileName

#aStarData, hdr = getdata(fileName, header=True)
#order = hdr['ORDERS'].replace(" ","")
#if "NIRSPEC" in hdr['INSTR']:
#    sourceName = hdr['TARGNAME'].replace(" ","")
##sourceName =fileName.split("_")[0]
#else:
#    sourceName = hdr['OBJECT'].replace(" ","")

#thismod = modelfileName.split('/')
#thismod = thismod[-1].split('.')
#suffix = '-dreamZPT-'
#if fsys: 
#    suffix += 'fsys-'
#suffix += thismod[0] + '-'

#xStar0 = aStarData[0,:]
#yStar0 = aStarData[1,:]
#errStar0 = aStarData[2,:]


#indir = sourceName + '-Results'
#outdir = indir

#def myHistogramPlot(x, n):
#    hist, bins = histogram(x, bins=n)
#    hist = 100.0 * hist / nansum(hist)
#    #    print nansum(hist)
#    width = bins[1] - bins[0]
#    center = (bins[:-1] + bins[1:])/2.0
#    bar(center, hist, align='center', width=width, color='gray')


infile = pklfileName
if os.path.exists(infile):
    pass
#else:
#    print 'Cannot find ', infile
#    suffix = suffix.replace('dreamZPT', 'dreamZS')
#    infile = indir + '/' + sourceName + '-' + order + suffix + 'Results.pkl.gz'
#    print 'Trying ', infile, ' instead'

if not os.path.exists(infile):
    print 'Error: cannot find ', infile
    sys.exit(1)

print infile
outfile = infile.replace('.pkl', 'vsini.txt')
print outfile

fp = zopen(infile, 'r')

db = load(fp)
fp.close()
lsf = db['FWHM'] * 1000.
alpha = db['alpha']
vsini = db['vsini']
cz = db['vr']
T = db['T']
logg = db['logg']
logp = db['ZP']
chisq = db['chisq']

burn = autoburn(chisq, 100)
if burn < 1000: burn = 1000
#print burn
chisq = chisq[burn:]
print 'Burn = ', burn
print chisq

# filter results
lsf = lsf[burn:]
alpha = alpha[burn:]
vsini = vsini[burn:]
print vsini
cz = cz[burn:]
T = T[burn:] 
logg = logg[burn:]

#np.savetxt(outfile,(lsf,alpha,cz,vsini,T,logg), fmt='%6.4f')
np.savetxt(outfile, vsini, fmt='%6.4f')

#f = open(outfile, 'w')


#print out the MC run to text

#f1 = ('%s %5.4f %5.4f %5.4f %5.4f %5.4f')
#f2 = ('%s %6.4f %6.4f %6.4f %6.4f %6.4f')
#f3 = ('%s %8.4f %8.4f %8.4f %8.4f %8.4f')
#f4 = ('%s %9.7f %9.7f %9.7f %9.7f %9.7f')
#f5 = ('%s %6.3f %6.3f %6.3f %6.3f %6.3f')
# sigma percentiles
# 0.1587
# 0.8414
#sap = lambda x, y: scoreatpercentile(x, y)
#print >>f,  'Parameter, median, 5%, 95% -s +s'
#print >>f,  '--------------------------------'
#print >>f,  f1 % ('    lsf (nm): ', median(lsf), sap(lsf, 5), sap(lsf, 95), \
#    median(lsf) - sap(lsf, 15.87), sap(lsf, 84.14) - median(lsf))
#print >>f,  f2 % ('      alpha : ', median(alpha), sap(alpha, 5), sap(alpha, 95), \
#    median(alpha) - sap(alpha, 15.87), sap(alpha, 84.14) - median(alpha))
#print >>f,  f5 % ('   cz (km/s): ', median(cz), sap(cz, 5), sap(cz, 95), \
#    median(cz) - sap(cz, 15.87), sap(cz, 84.14) - median(cz))
#print >>f,  f5 % ('vsini (km/s): ', median(vsini), sap(vsini, 5), sap(vsini, 95), \
#    median(vsini) - sap(vsini, 15.87), sap(vsini, 84.14) - median(vsini))
#print >>f,  f5 % ('T surface(K): ', median(T), sap(T, 5), sap(T, 95), \
#    median(T) - sap(T, 15.87), sap(T, 84.14) - median(T))
#print >>f,  f5 % ('      log(g): ', median(logg), sap(logg, 5), sap(logg, 95), \
#                  median(logg) - sap(logg, 15.87), sap(logg, 84.14) - median(logg))
#'''
#print >>f,  f4 % ('   l0 (mu m): ', median(l0), sap(l0, 5), sap(l0, 95), \
#    median(l0) - sap(l0, 15.87), sap(l0, 84.14) - median(l0))
#print >>f,  f4 % ('   l1 (mu m): ', median(l1), sap(l1, 5), sap(l1, 95), \
#    median(l1) - sap(l1, 15.87), sap(l1, 84.14) - median(l1))
#print >>f,  f4 % ('   l2 (mu m): ', median(l2), sap(l2, 5), sap(l2, 95), \
#    median(l2) - sap(l2, 15.87), sap(l2, 84.14) - median(l2))
#'''
#print >>f, 'Minimum chisq = ', amin(chisq)
#print >>f, '   Red. chisq = ', amin(chisq) / db['dof']
##print >>f, '          rms = ', rms
#print >>f, 'Median resid. = ', nanmedian(sqrt(res**2))
#print >>f, 'Median uncer. = ', nanmedian(sqrt(dyStar**2))
#print >>f, 'System uncer. = ', sysunc
#f.close()


#f1 = ('%s %5.4f %5.4f')
#f2 = ('%s %6.4f %6.4f')
#f3 = ('%s %8.4f %8.4f')
#f4 = ('%s %9.7f %9.7f')
#f5 = ('%s %6.3f %6.3f')

#print >>f,  'Parameter, mean, stdev'
#print >>f,  '--------------------------------'
#print >>f,  f1 % ('    lsf (nm): ', mean(lsf), std(lsf))
#print >>f,  f2 % ('      alpha : ', mean(alpha), std(alpha))
#print >>f,  f5 % ('   Vr (km/s): ', mean(cz), std(cz))
#print >>f,  f5 % ('vsini (km/s): ', mean(vsini), std(vsini))
#print >>f,  f5 % ('T surface(K): ', mean(T), std(T))
#print >>f,  f5 % ('      log(g): ', mean(logg), std(logg))
#'''
#print >>f,  f4 % ('   l0 (mu m): ', mean(l0), std(l0))
#print >>f,  f4 % ('   l1 (mu m): ', mean(l1), std(l1))
#print >>f,  f4 % ('   l2 (mu m): ', mean(l2), std(l2))
#'''
#print >>f, 'Minimum chisq = ', amin(chisq)
#print >>f, '   Red. chisq = ', amin(chisq) / db['dof']
#print >>f, '          rms = ', rms
#f.close()

