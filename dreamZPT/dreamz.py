'''
Jack Gallimore, Bucknell University, 2015
'''
from pylab import *
from random import sample as sampleWithoutReplacement


x0 = ones(2) * 9.0
x1 = ones(2) * 4.0
s0 = 0.1
s1 = 0.1

def twoPeaklogP(p):
    arg0 = -dot(p - x0, p - x0) / s0 / s0
    arg1 = -dot(p - x1, p - x1) / s1 / s1
    return logaddexp(arg0, arg1)

def priorP(p):
    if p[0] < 0: return 0
    if p[0] > 20: return 0
    if p[1] < 0: return 0
    if p[1] > 20: return 0
    return 1.0

def initializeDREAM(logPfunc, nChains, ndim, ninit):
    Z = uniform(size=(ninit, ndim)) * 20.
    ZP = zeros(ninit)
    for i in range(ninit):
        ZP[i] = logPfunc(Z[i])
    X = uniform(size=(nChains, ndim)) * 20.
    XP = zeros(nChains)
    for i in range(nChains):
        XP[i] = logPfunc(X[i])
    return X, XP, Z, ZP


def dreamzStep(func, X, XP, Z, ZP, gamma=0.1, eps=1.e-6, e=0.05):
    nChains, ndim = shape(X)
    nHist = len(ZP)
    idx = arange(nHist)
    for i in range(nChains):
        r1, r2 = sampleWithoutReplacement(idx, 2)
        delta = Z[r2] - Z[r1]
        x = X[i] + (1.0+normal()*e)*gamma*delta + normal()*eps
        xp = func(x)
        alpha = exp(xp - XP[i])
        if uniform() < alpha:
            X[i] = x
            XP[i] = xp
    return X, XP
        

nChains = 5
niter = 10000
thin = 10

X, XP, Z, ZP = initializeDREAM(twoPeaklogP, nChains, 2, 100)

for i in range(niter):
    for j in range(thin):
        X, XP = dreamzStep(twoPeaklogP, X, XP, Z, ZP)
    Z = append(Z, X, axis=0)
    ZP = append(ZP, XP)

