import os.path

trange = array([12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, \
                17, 17.5, 18, 18.5, 19, 19.5, 20, 21, 22, 23, 24, 25, 26])
tList = []
for t in trange:
    if t - floor(t) == 0.5:
        thing = 'lte%05.1f' % (t,)
    else:
        thing = 'lte%03i' % (t,)
    tList += [thing]

gRange = arange(3.5, 6.0, 0.5)
gList = []
for g in gRange:
    gstr = '%3.1f' % (g,)
    gList += [gstr]

# check if files exist
for t in tList:
    for g in gList:
        infile = t + '-' + g + '-0.0a+0.0.BT-Settl.spec.7.bz2'
        if not(os.path.exists(infile)):
            print infile

'''
Missing files:
lte014.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2
lte016.5-3.5-0.0a+0.0.BT-Settl.spec.7.bz2

Correspond to cube indices
(:, 5, 1)
(:, 9, 0)
'''

# ok, let's assemble a cube
# load model atmosphere
def expDtofloat(s): # deal with FORTRAN "D"s *shiver*
    return float(s.replace('D', 'E'))

def loadspectrum(infile):
    x, y = loadtxt(infile, unpack=True, \
                   usecols=(0,1),\
                   converters={0:expDtofloat, 1:expDtofloat})
    x = x / 10000.
    y = 10.**(y-12.)
    
    idx = logical_and(x >= 2.0, x <= 2.4)
    x = x[idx]
    y = y[idx]

    return(x, y)

x, y = loadspectrum(infile)
nx = len(x)

ny = len(trange)
nz = len(gRange)

modelCube = zeros((nx, ny, nz))

for i in range(len(tList)):
    t = tList[i]
    for j in range(len(gList)):
        g = gList[j]
        print t, g
        infile = t + '-' + g + '-0.0a+0.0.BT-Settl.spec.7.bz2'
        if os.path.exists(infile):
            x, y = loadspectrum(infile)
        else:
            y = ones(40001) * -999.
        modelCube[:,i, j] = y

# patch the cube for missing models
patch0 = modelCube[:,4,1]
patch1 = modelCube[:,6,1]
#patch2 = modelCube[:,5,0]
#patch3 = modelCube[:,5,2]

#modelCube[:,5,1] = 0.25 * (patch0 + patch1 + patch2 + patch3)
modelCube[:,5,1] = 0.5 * (patch0 + patch1)

patch0 = modelCube[:,8,0]
patch1 = modelCube[:,10,0]
modelCube[:,9,0] = 0.5 * (patch0 + patch1)


db = {}
db['wave'] = x
db['T'] = trange
db['g'] = gRange
db['models'] = modelCube

from cPickle import dump
from gzip import open as zopen

outfile = 'btsettl.pkl.gz'
f = zopen(outfile, 'w')
dump(db, f)
f.close()


# need to repair missing files
