'''
Jack Gallimore, Bucknell University, 2015
'''
from scipy import *
from numpy import *
from numpy.random import normal
from numpy.random import uniform
#from numpy.random import rand as uniform
from random import sample as sampleWithoutReplacement

def sampleIndexWithoutReplacement(N, k, mask = [-1]):
    result = zeros(k,dtype=int) - 1
    for i in range(k):
        j = int(N*uniform())
        while j in result or j in mask:
            j = int(N*uniform())
        result[i] = j
    return result


'''
Perform differential evolution (DREAM - style).
'''
def gwStep(func, priorFunc, \
           X, XP, \
           scale=2, 
           fitVector=None, args=None):
    antiScale = 1./scale
    nChains, ndim = shape(X)
    nHist = len(XP)
    XPNEW = copy(XP)
    XNEW = copy(X)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    for i in range(nChains):
        ''' 
        Draw another index at random
        '''
        r = sampleIndexWithoutReplacement(nHist, 1, mask=[i])[0]
        #Z = sqrt(uniform()*(scale - antiScale) + antiScale)
        Z = ((scale - 1.)*uniform()+1)**2/scale
        '''
        GW step: make a trial, new solution.
        '''
        xk = copy(X[i])
        xj = X[r]
        xk[parIndex] = xj[parIndex] + Z*(xk[parIndex]-xj[parIndex])
        pf = priorFunc(xk)
        if ~isneginf(pf):
            if args == None:
                xp = func(xk) + pf
            else:
                xp = func(xk, *args) + pf
        else:
            xp = 0. # doesn't matter

        '''
        Acceptance ratio setp. 
        '''
        alpha = xp - XP[i] +pf +(ndim-1.)*log(Z)
        dice = log(uniform())
        if dice < alpha:
            XNEW[i] = xk
            XPNEW[i] = xp
            success += [1]
        else:
            success += [0]
    return XNEW, XPNEW, success



'''
Perform differential evolution (DREAMZ-style).
'''
def gwStepZ(func, priorFunc, \
           X, XP, Z, \
           scale=2, 
           fitVector=None, args=None):
    antiScale = 1./scale
    nChains, ndim = shape(X)
    nHist = shape(Z)[0]
    XPNEW = copy(XP)
    XNEW = copy(X)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    for i in range(nChains):
        ''' 
        Draw another index at random
        '''
        r = sampleIndexWithoutReplacement(nHist, 1)[0]
        z = ((scale - 1.)*uniform()+1)**2/scale
        '''
        GW step: make a trial, new solution.
        '''
        xk = copy(X[i])
        xj = Z[r]
        xk[parIndex] = xj[parIndex] + z*(xk[parIndex]-xj[parIndex])
        pf = priorFunc(xk)
        if ~isneginf(pf):
            if args == None:
                xp = func(xk) + pf
            else:
                xp = func(xk, *args) + pf
        else:
            xp = 0. # doesn't matter

        '''
        Acceptance ratio setp. 
        '''
        alpha = xp - XP[i] +pf +(ndim-1.)*log(z)
        dice = log(uniform())
        if dice < alpha:
            XNEW[i] = xk
            XPNEW[i] = xp
            success += [1]
        else:
            success += [0]
    return XNEW, XPNEW, success


