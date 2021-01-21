#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as mp
mp.ioff()
from pylab import *
mpl.rc('font',family='Times New Roman')
from pylab import *
from glob import glob
from scipy.stats import scoreatpercentile
from astrolib.FindHelio import findCorr
from dreamZPT.diagnostics import autoburn
import argparse
from cPickle import load
from gzip import open as zopen
from pyfits import getdata
import matplotlib

parser = argparse.ArgumentParser(description='Plot a fit to a Brown Dwarf Spectrum')
parser.add_argument('sourceName', help='Name of Brown Dwarf (prefix to Pickle filename)')
parser.add_argument('--fsys', help='Whether or not systemic uncertainties were included. Toggles True (default=False)', default=False, action='store_true')
parser.add_argument('--modelfile', help='Name of model grid (default: ./CIFIST2011/btsettl-4.pkl.gz', default='./CIFIST2011/btsettl-4.pkl.gz')

args = parser.parse_args()
sourceName = args.sourceName
fsys = args.fsys
modelfileName = args.modelfile
thismod = modelfileName.split('/')
thismod = thismod[-1].split('.')
suffix = '-dreamZPT-'
if fsys: 
    suffix += 'fsys-'
suffix += thismod[0] + '-'

indir = './' + sourceName + '-Results/'
wildCard = indir + '*-??-dreamZS-'
if fsys: wildCard += 'fsys-'
wildCard += thismod[0] + '-Results.pkl.gz'

fl = sorted(glob(wildCard))

print wildCard
print len(fl)

if len(fl) == 0:
    wildCard = indir + '*-??-dreamZPT-'
    if fsys: wildCard += 'fsys-'
    wildCard += thismod[0] + '-Results.pkl.gz'
    fl = sorted(glob(wildCard))

print fl

vrList = []
vsiniList = []
chiList = []
fwhmList = []
TList = []
gList = []
alphaList = []

for f in fl:
    fp = zopen(f, 'r')
    db = load(fp)
    fp.close()
    b = db['burn']
    if b < 1000: b = 1000
    vrList += [db['vr'][b:]]
    vsiniList += [db['vsini'][b:]]
    chiList += [db['chisq'][b:]/db['dof']]
    fwhmList += [db['FWHM'][b:]*1000.]
    TList += [db['T'][b:]]
    gList += [db['logg'][b:]]
    alphaList += [db['alpha'][b:]]

clf()
boxplot(vrList, whis='range')
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('RV (km/s)')
xlabel('Order')

savefig(indir + sourceName + '-vr-boxplot.pdf')

clf()
boxplot(vsiniList, whis='range')
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('$v$ sin($i$) (km/s)')
xlabel('Order')

savefig(indir + sourceName + '-vsini-boxplot.pdf')

clf()
boxplot(chiList, whis='range')
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('$\\chi^2$ (reduced)')
xlabel('Order')

savefig(indir + sourceName + '-chi2-boxplot.pdf')

clf()
boxplot(fwhmList, whis='range')
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('LSF FWHM (nm)')
xlabel('Order')

savefig(indir + sourceName + '-fwhm-boxplot.pdf')

clf()
boxplot(TList, whis='range')
print TList
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('$T_{eff}$ (K)')
xlabel('Order')

savefig(indir + sourceName + '-T-boxplot.pdf')

clf()
boxplot(gList, whis='range')
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('log(g)')
xlabel('Order')

savefig(indir + sourceName + '-logg-boxplot.pdf')

clf()
boxplot(alphaList, whis='range')
xticks(arange(7)+1, ('32', '33', '34', '35', '36', '37', '38'))

ylabel('$\\alpha$ (Telluric Scaling)')
xlabel('Order')

savefig(indir + sourceName + '-alpha-boxplot.pdf')

