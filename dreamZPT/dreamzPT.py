'''
Jack Gallimore, Bucknell University, 2015
'''
from pylab import *
from random import sample as sampleWithoutReplacement


x0 = ones(2) * 9.0
x1 = ones(2) * 4.0
s0 = 0.1
s1 = 0.1

# log probability function - after Maggie
def twoPeaklogP(p):
    arg0 = -dot(p - x0, p - x0) / s0 / s0
    arg1 = -dot(p - x1, p - x1) / s1 / s1
    return logaddexp(arg0, arg1)

# Prior probabilities on parameters -- not used yet
def priorP(p):
    if p[0] < 0: return 0
    if p[0] > 20: return 0
    if p[1] < 0: return 0
    if p[1] > 20: return 0
    return 1.0

'''
Initialize Z and X arrays with random values between 0 and 20. Each
chain in X and Z are assigned a temperature. X[0] is the low
temperature (T = 1) chain. X[1] has a temperature of 2, X[2] has a
temperature of 4, etc. To get higher temperatures in the run, raise
the number of chains.

XP, ZP store probabilities
XT, ZT store temperatures
'''

def initializeDREAMPT(logPfunc, nChains, ndim, ninit):
    Z = uniform(size=(ninit, ndim)) * 20.
    ZP = zeros(ninit)
    ZT = 2.0**mod(arange(ninit), nChains)
    for i in range(ninit):
        ZP[i] = logPfunc(Z[i])
    X = uniform(size=(nChains, ndim)) * 20.
    XT = 2.0**arange(nChains)
    XP = zeros(nChains)
    for i in range(nChains):
        XP[i] = logPfunc(X[i])
    return X, XP, XT, Z, ZP, ZT

'''
Perform differential evolution (DREAMZ). Parallel tempering
modification: only allow single temperature evolution. In other words,
a given X[i] only is modified by Z values at the same temperature.
'''
def dreamzPTStep(func, X, XP, XT, Z, ZP, ZT, gamma=1.0, eps=1.e-6, e=0.05):
    nChains, ndim = shape(X)
    nHist = len(ZP)
    #idx = arange(nHist)
    for i in range(nChains):
        '''
        Make a list of Z index values at the same temperature as X[i].
        '''
        idx = arange(nHist)[ZT==XT[i]]
        '''
        Sample two indices at random.
        '''
        r1, r2 = sampleWithoutReplacement(idx, 2) # select random
        '''
        Create an evolution vector from those samples.
        '''
        delta = Z[r2] - Z[r1]
        '''
        Differential evolution step: make a trial, new solution.
        '''
        x = X[i] + (1.0+normal()*e)*gamma*delta + normal()*eps
        xp = func(x)
        '''
        Acceptance ratio setp. Notice the inclusion of temperature.
        '''
        alpha = exp((xp - XP[i])/XT[i])
        if uniform() < alpha:
            X[i] = x
            XP[i] = xp
    return X, XP

'''
Swap step: this function allows low temperature and high temperature
solutions to talk to each other. Take a random parameter vector in X
and try swapping it with the next higher temperature parameter vector
in X. The acceptance ratio is based on Earl and Deem (2008). 
'''
def dreamPTSwap(func, X, XP, XT):
    nChains, ndim = shape(X)
    XNEW = copy(X)
    XPNEW = copy(XP)
    i = sampleWithoutReplacement(range(1, nChains), 1)[0]
    j = i-1
    betaLow = 1.0 / XT[j]
    betaHigh = 1.0 / XT[i]
    alpha = exp((XP[i] - XP[j])*(betaHigh - betaLow))
    if uniform() < alpha:
        XNEW[i] = X[j]
        XNEW[j] = X[i]
        XPNEW[i] = XP[j]
        XPNEW[j] = XP[i]
    return XNEW, XPNEW

nChains = 5
niter = 10000
thin = 10

X, XP, XT, Z, ZP, ZT= initializeDREAMPT(twoPeaklogP, nChains, 2, 100)

for i in range(niter):
    # Perform differential evolution for "thin" iterations
    for j in range(thin):
        X, XP = dreamzPTStep(twoPeaklogP, X, XP, XT, Z, ZP, ZT)
    # At the end of thinning, try a random swap
    X, XP = dreamPTSwap(twoPeaklogP, X, XP, XT)
    # Update the history (Z) arrays
    Z = append(Z, X, axis=0)
    ZP = append(ZP, XP)
    ZT = append(ZT, XT)

# write out the results

nhist, ndim = shape(Z)
ZZ = zeros((nhist, ndim+2))
ZZ[:,0:ndim] = Z
ZZ[:,ndim] = ZP
ZZ[:,ndim+1] = ZT
savetxt('ptTest.txt', ZZ)

# plot the low temperature results (in blue) on top of the highest temperature
# results (in red)
fg = figure(num=1, figsize=(8,8))
clf()
plot(Z[ZT==XT[-1], 0], Z[ZT==XT[-1], 1], 'r,')
plot(Z[ZT==XT[0], 0], Z[ZT==XT[0], 1], 'b,')
xlim(0,10)
ylim(0,10)
xlabel('Parameter 1')
ylabel('Parameter 2')
savefig('ptTest.png')

