'''
Jack Gallimore, Bucknell University, 2015
'''
from scipy import *
from numpy import *
from numpy.random import normal
from numpy.random import uniform
#from numpy.random import rand as uniform
from random import sample as sampleWithoutReplacement
from multiprocessing.dummy import Pool # multithreading
#from multiprocessing import Pool # multiprocessing
from functools import partial

def sampleIndexWithoutReplacement(N, k, mask = [-1]):
    result = zeros(k,dtype=int) - 1
    for i in range(k):
        j = int(N*uniform())
        while j in result or j in mask:
            j = int(N*uniform())
        result[i] = j
    return result
    

'''
Perform differential evolution (DREAMZ). Parallel tempering
modification: only allow single temperature evolution. In other words,
a given X[i] only is modified by Z values at the same temperature.

Note: assumes 1 temperature per chain
'''
def dreamzPTStep(func, priorFunc, \
                     X, XP, XT, Z, ZP, ZT, \
                     ncr=0, navg=1, \
                     gamma=1.0, eps=1.e-6, e=0.05, \
                     fitVector=None, args=None, plow=None, phigh=None):
    nChains, ndim = shape(X)
    nHist = len(ZP)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    #idx = arange(nHist)
    success = list(zeros(nChains))
    idx = arange(nHist)
    #mask = zeros(nHist, dtype=bool)
    nn = nHist / nChains
    XNEW = X * 0.0 
    XPNEW = XP * 0.0
    '''
    sfunc = partial(getLogProb, func, priorFunc, \
                    args)
    '''
    def evolve(i):
        '''
        Make a list of Z index values at the same temperature as X[i].
        '''
        thisTemp = XT[i]
        '''
        Sample indices at random. Take into account averaging.
        '''
        r = sampleIndexWithoutReplacement(nn,2*navg)*nChains + i

        Z1 = sum(Z[r[0:navg],:], axis=0)
        Z2 = sum(Z[r[navg:],:], axis=0)
        '''
        Create an evolution vector from those samples.
        '''
        delta = Z2 - Z1
        '''
        Differential evolution step: make a trial, new solution.
        '''
        x = copy(X[i])
        pp = copy(parIndex)
        if ncr > 0 and ncr < ndim: # deal with limited number of crossovers
            pp = sampleWithoutReplacement(parIndex, ncr)
        x[pp] += (1.0+normal(size=len(pp))*e)*gamma*delta[pp] + normal(size=len(pp))*eps
        ''' 
        Try reflecting
        '''
        if plow != None and phigh != None:
            idx = (x > phigh) 
            x[idx] = 2*phigh[idx] - x[idx]
            idx = (x < plow)
            x[idx] = 2*plow[idx] - x[idx]
        pf = priorFunc(x)
        if ~isneginf(pf):
            if args == None:
                xp = func(x)
            else:
                xp = func(x, *args)
        else:
            xp = 0. # doesn't matter
        xp = xp + pf
        alpha = (xp - XP[i])/XT[i]
        dice = log(uniform())
        successi = dice < alpha
        if successi:
            X[i] = x
            XP[i] = xp
        success[i] = int(successi)
        return
    
    pool = Pool(processes=3)
    pool.map(evolve, range(nChains))
    pool.close()
    pool.join()

    #success = list(success.astype(int))
    return X, XP, success

'''
Swap step: this function allows low temperature and high temperature
solutions to talk to each other. Take a random parameter vector in X
and try swapping it with the next higher temperature parameter vector
in X. The acceptance ratio is based on Earl and Deem (2008). 
'''
def dreamzPTSwap(func, priorFunc, X, XP, XT, args=None):
    TU = unique(XT)
    ntemps = len(TU)
    nChains, ndim = shape(X)
    #XNEW = copy(X)
    #XPNEW = copy(XP)
    for i in range(len(TU)-1,0,-1):
        Tlow = TU[i-1]
        Thigh = TU[i]
        betaLow = 1.0 / Tlow
        betaHigh = 1.0 / Thigh
        #tfactor = (Tlow - Thigh)/(Tlow*Thigh)

    # pick chains with matching temperatures to swap
        idx = where(XT == Tlow)[0]
        j = sampleWithoutReplacement(idx, 1)[0]
        Plow = XP[j]
        idx = where(XT == Thigh)[0]
        i = sampleWithoutReplacement(idx, 1)[0]
        Phigh = XP[i]

        '''
        pfi = priorFunc(X[i])
        pfj = priorFunc(X[j])
        '''

        alpha = (Phigh - Plow)*(betaLow - betaHigh)
        if log(uniform()) < alpha:
            xx = copy(X[j])
            xxp = copy(XP[j])
            X[j] = X[i]
            XP[j] = XP[i]
            X[i] = xx
            XP[i] = xxp
    return X, XP

def getLogProb(func, priorFunc, args, p):
    xp = 0.
    pf = priorFunc(p)
    if ~isneginf(pf):
        if args == None:
            xp = func(p)
        else:
            xp = func(p, *args)
    return xp + pf



'''
Perform differential evolution (DREAMZ).
Added reflections at boundaries if plow and phigh are defined.
Plow = vector of lower boundary values, one for each parameter.
Phigh = ditto for upper boundary values.
'''
def dreamzStep(func, priorFunc, \
                   X, XP, Z, ZP, \
                   ncr=0, navg=1, \
                   gamma=1.0, eps=1.e-6, e=0.05, \
                   fitVector=None, args=None, plow=None, phigh=None):

    nChains, ndim = shape(X)
    nHist = len(ZP)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    #idx = arange(nHist)
    XNEW = X * 0.0
    sfunc = partial(getLogProb, func, priorFunc, \
                    args)
    for i in range(nChains):
        '''
        Sample two indices at random.
        '''
        #r = sampleWithoutReplacement(idx, 2*navg) # select random
        r = sampleIndexWithoutReplacement(nHist, 2*navg)
        Z1 = sum(Z[r[0:navg],:], axis=0)
        Z2 = sum(Z[r[navg:],:], axis=0)
        '''
        Create an evolution vector from those samples.
        '''
        delta = Z2 - Z1
        '''
        Differential evolution step: make a trial, new solution.
        '''
        x = copy(X[i])
        pp = copy(parIndex)
        if ncr > 0 and ncr < ndim: # deal with limited number of crossovers
            pp = sampleWithoutReplacement(parIndex, ncr)
        x[pp] += (1.0+normal(size=len(pp))*e)*gamma*delta[pp] + \
                 normal(size=len(pp))*eps
        ''' 
        Try reflecting
        '''
        if plow != None and phigh != None:
            idx = (x > phigh) 
            x[idx] = 2*phigh[idx] - x[idx]
            idx = (x < plow)
            x[idx] = 2*plow[idx] - x[idx]
        XNEW[i] = x

    #pool = Pool(processes=2)   
    pool = Pool()
    XPNEW = array(pool.map(sfunc, XNEW))
    pool.close()
    pool.join()
    #print XPNEW
    alpha = XPNEW - XP
    dice = log(uniform(size=nChains))
    success = dice < alpha
    XP[success] = XPNEW[success]
    X[success,:] = XNEW[success,:]
    success = list(success.astype(int))
    return X, XP, success

'''
Perform differential evolution (DREAM).
'''
def dreamStep(func, priorFunc, \
                   X, XP, \
                   ncr=0, navg=1, \
                   gamma=1.0, eps=1.e-6, e=0.05, \
                   fitVector=None, args=None, violateDB=False):

    nChains, ndim = shape(X)
    nHist = len(XP)
    XPNEW = copy(XP)
    XNEW = copy(X)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    #idx = arange(nHist)
    #mask = ones(nHist, dtype=bool)
    for i in range(nChains):
        '''
        Sample two indices at random.
        '''
        #if not violateDB: # violateDB may be desirable during burn-in
            #mask[i] = False # ensure current chain is not included in selection
        if violateDB:
            r = sampleIndexWithoutReplacement(nHist, 2*navg)
        else:
            r = sampleIndexWithoutReplacement(nHist, 2*navg, mask=[i])
        #r = sampleWithoutReplacement(idx[mask], 2*navg) # select random
        #mask[i] = True # restore current chain for future iterations
        Z1 = sum(X[r[0:navg],:], axis=0)
        Z2 = sum(X[r[navg:],:], axis=0)
        '''
        Create an evolution vector from those samples.
        '''
        delta = Z2 - Z1
        '''
        Differential evolution step: make a trial, new solution.
        '''
        x = copy(X[i])
        pp = copy(parIndex)
        if ncr > 0 and ncr < ndim: # deal with limited number of crossovers
            pp = sampleWithoutReplacement(parIndex, ncr)
        x[pp] += (1.0+normal(size=len(pp))*e)*gamma*delta[pp] + \
                 normal(size=len(pp))*eps
        pf = priorFunc(x)
        if ~isneginf(pf):
            if args == None:
                xp = func(x)
            else:
                xp = func(x, *args)
        else:
            xp = 0. # doesn't matter

        '''
        Acceptance ratio setp. 
        '''
        alpha = xp + pf - XP[i]
        dice = log(uniform())
        if dice < alpha:
            XNEW[i] = x
            XPNEW[i] = xp + pf
            success += [1]
        else:
            success += [0]
    return XNEW, XPNEW, success

'''
Perform differential evolution (standard -- violates detailed balance.
'''
def DEStep(func, priorFunc, \
                   X, XP, \
                   ncr=0, \
                   gamma=1.0, eps=1.e-6, e=0.05, \
                   fitVector=None, args=None):

    nChains, ndim = shape(X)
    nHist = len(XP)
    XPNEW = copy(XP)
    XNEW = copy(X)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    idx = arange(nHist)
    mask = ones(nHist, dtype=bool)
    for i in range(nChains):
        '''
        Sample three indices at random.
        '''
        mask[i] = False # ensure current chain is not included in selection
        r = sampleWithoutReplacement(idx[mask], 3) # select random
        mask[i] = True # restore current chain for future iterations
        Z0 = X[r[0]]
        Z1 = X[r[1]]
        Z2 = X[r[2]]
        '''
        Create an evolution vector from those samples.
        '''
        delta = Z2 - Z1
        '''
        Differential evolution step: make a trial, new solution.
        '''
        x = copy(X[i])
        pp = copy(parIndex)
        if ncr > 0 and ncr < ndim: # deal with limited number of crossovers
            pp = sampleWithoutReplacement(parIndex, ncr)
        x[pp] = Z0[pp] + (1.0+normal(size=len(pp))*e)*gamma*delta[pp] + \
                 normal(size=len(pp))*eps
        pf = priorFunc(x)
        if ~isneginf(pf):
            if args == None:
                xp = func(x)
            else:
                xp = func(x, *args)
        else:
            xp = 0. # doesn't matter

        '''
        Acceptance ratio setp. 
        '''
        alpha = xp + pf - XP[i]
        dice = log(uniform())
        if dice < alpha:
            XNEW[i] = x
            XPNEW[i] = xp + pf
            success += [1]
        else:
            success += [0]
    return XNEW, XPNEW, success

'''
Perform differential evolution (DREAMZ).
'''
def dreamSwap(func, priorFunc, \
              X, XP, \
              eps=1.e-6, 
              fitVector=None, args=None):

    nChains, ndim = shape(X)
    nHist = len(XP)
    XPNEW = copy(XP)
    XNEW = copy(X)
    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    idx = arange(nHist)
    mask = ones(nHist, dtype=bool)
    for i in range(nChains):
        '''
        Sample another index at random.
        '''
        mask[i] = False # ensure current chain is not included in selection
        r = sampleWithoutReplacement(idx[mask], 1) # select random
        mask[i] = True # restore current chain for future iterations
        Z1 = X[r[0],:]

        x = copy(X[i])
        pp = copy(parIndex)
        x[pp] = Z1[pp] + normal(size=len(pp))*eps
        pf = priorFunc(x)
        if ~isneginf(pf):
            if args == None:
                xp = func(x)
            else:
                xp = func(x, *args)
        else:
            xp = 0. # doesn't matter

        '''
        Acceptance ratio setp. 
        '''
        alpha = xp + pf - XP[i]
        dice = log(uniform())
        if dice < alpha:
            XNEW[i] = x
            XPNEW[i] = xp + pf
    return XNEW, XPNEW

from numpy.linalg import norm as vnorm

def vecProj(z, x): # project vector z onto x. Returns projected vector
    nrm = vnorm(x,2)
    unitVector = x / nrm
    return dot(z,x)*unitVector / nrm

def snookerStep(func, priorFunc, \
                    X, XP, Z, \
                    fitVector=None, args=None):
    nChains, ndim = shape(X)
    M, nparms = shape(Z)
    XNEW = copy(X)
    XPNEW = copy(XP)
    bias = copy(XPNEW)

    if fitVector == None:
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)

    d = len(parIndex)
    #idz = arange(M)
    success = []
    sfunc = partial(getLogProb, func, priorFunc, \
                    args)
    for i in range(nChains):
        x = copy(X[i])
        xp = XP[i]
        # now pick random Z's for snooker
        idx = sampleIndexWithoutReplacement(M, 3)
        z = Z[idx[0]]
        xMz = x - z
        if abs(sum(xMz)) > 0:
            # project onto x-z
            zr1 = vecProj(Z[idx[1]], xMz)
            zr2 = vecProj(Z[idx[2]], xMz)
            gamma = uniform()*(2.2 - 1.2) + 1.2
            diff = gamma * (zr1 - zr2)
            #x[parIndex] = x[parIndex] + diff[parIndex]
            XNEW[i, parIndex] = x[parIndex] + diff[parIndex]
            bias[i] = pow(vnorm((x-z), 2) / \
                       vnorm((X[i] - z), 2), d-1)
            

    pool = Pool(processes=2)   
    XPNEW = array(pool.map(sfunc, XNEW))
    pool.close()
    pool.join()
    alpha = XPNEW - XP + log(bias)
    dice = log(uniform(size=nChains))
    success = dice < alpha
    XP[success] = XPNEW[success]
    X[success,:] = XNEW[success,:]
    success = list(success.astype(int))
    return X, XP, success

