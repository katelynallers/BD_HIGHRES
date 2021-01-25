from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d, convolve1d, \
    uniform_filter1d
from cPickle import load, dump
from gzip import open as zopen
c = 299792.458 # km/s
OVRSAMPL = 10.

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

fp = zopen('btsettl.pkl.gz')
db = load(fp)
fp.close()

w = db['wave']
w0 = log(amin(w))
w1 = log(amax(w))
mwlVgrid = exp(linspace(w0, w1, len(w)*OVRSAMPL))
vpix = c * median(diff(mwlVgrid)/mwlVgrid[1:])

m = db['models']
T = db['T']
g = db['g']

vmax = 60.
vmin = 0.
dv = 5.

vrange = arange(vmin,vmax+dv, dv)
nv = len(vrange)

mOut = zeros((60001, 23, 5, nv))

for iv in range(nv):
    vsini = vrange[iv]
    for it in range(len(T)):
        for ig in range(len(g)):
            model = m[:,it,ig]
            # convolve with vsini
            mflVgrid = InterpolatedUnivariateSpline(w, model, k=3)(mwlVgrid)
            velgrid, lsfRot = lsf_rotate(vpix, vsini)
            if vsini > vpix:
                mflVgridC = convolve1d(mflVgrid, lsfRot, mode='nearest')
            else:
                mflVgridC = copy(mflVgrid)
            cmodel = InterpolatedUnivariateSpline(mwlVgrid, mflVgridC, k=3)(w)
            mOut[:,it,ig,iv] = cmodel

db['models'] = mOut
db['vsini'] = vrange

fp = zopen('btsettl-vsini.pkl.gz', 'w')
dump(db, fp)
fp.close()

            
