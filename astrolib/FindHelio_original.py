#!/usr/bin/env python
import numpy as np
from numpy import sin, cos
from pyfits import getheader
#from astrolib.helcorr import helcorr as findHelio
from astrolib.baryvel import baryvel as findVelo
import sys
radFrac = 180.0/ np.pi#change to radians


#infile = sys.argv[1] + '.fits'#name of file from terminal
#m = getheader(infile)#read
#print m
'''
#print l 
#print m[25]
#print m[26]
RA = (m['RA'])/radFrac
MRA = RA/15
DEC = (m['DEC'])/radFrac
MJD = m['MJDOBS']
JD = MJD + 2400000.5#convert to julian date

#barrCorr, hjd, velo = findHelio(155.4681, 19.8207, 4205, MRA, DEC, JD, debug=False)

#print m
#print RA
#print DEC
#print barrCorr
#print hjd
#print m
#print JD

hVelo, bVelo = findVelo(JD, 2000)
print hVelo
print bVelo

vhelio = bVelo[0]*cos(DEC)*cos(RA) + bVelo[1]*cos(DEC)*sin(RA) + bVelo[2]*sin(DEC)
print vhelio
'''

def findCorr(infile):
    m = getheader(infile)#read
    RA = (m['RA'])/radFrac
    DEC = (m['DEC'])/radFrac
    try:
        MJD = m['MJD-OBS']
    except:
        pass
    try:
        MJD = m['MJD_OBS']
    except:
        pass
    try:
        MJD = m['MJDOBS']
    except:
        pass

    JD = MJD + 2400000.
    hVelo, bVelo = findVelo(JD, 2000)
    vhelio = bVelo[0]*cos(DEC)*cos(RA) + bVelo[1]*cos(DEC)*sin(RA) + bVelo[2]*sin(DEC)
    return vhelio


