'''
Jack Gallimore, Bucknell University, 2015
'''
from numpy import shape, arange, copy, zeros, isneginf, sum
from numpy import log, unique, where, ones, dot, inf, isfinite
from numpy import cumsum, exp, sqrt
from numpy.random import normal, uniform, randint
from random import sample as sampleWithoutReplacement

def sampleIndexWithoutReplacement(N, k, mask = [-1]):
    result = zeros(k,dtype=int) - 1
    for i in range(k):
        j = int(N*uniform())
        while j in result or j in mask:
            j = int(N*uniform())
        result[i] = j
    return result
    
def ncrDraw(ncrScl, ncrMax):
    selectFrom = arange(ncrMax) + 1
    arg = (1.0*selectFrom) / (1.0 * ncrScl)
    probs = exp(-arg*arg)
    cdf = cumsum(probs) / sum(probs)
    dice = uniform()
    ncr = selectFrom[where(dice < cdf)[0][0]]
    return ncr

'''
Perform differential evolution (DREAMZ). Parallel tempering
modification: only allow single temperature evolution. In other words,
a given X[i] only is modified by Z values at the same temperature.

Note: assumes 1 temperature per chain
'''
def dreamzPTStep(func, priorFunc, \
                     X, XP, XT, Z, ZP, ZT, \
                     ncr=0, navg=1, \
                     eps=1.e-6, e=0.05, \
                     fitVector=None, args=None, plow=None, phigh=None, thin=10):
    nChains, ndim = shape(X)
    nHist = len(ZP)
    if all(fitVector == None):
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    #idx = arange(nHist)
    success = []
    idx = arange(nHist)
    #mask = zeros(nHist, dtype=bool)
    nn = nHist / nChains
    if any(fitVector != None):
        nfit = len(fitVector)
    else:
        nfit = ndim    
    
    for j in range(thin):
        for i in range(nChains):
            nav = randint(1,navg+1)
            if ncr > 0: 
                #nc = np.random.randint(1,ncr+1)
                nc = ncrDraw(ncr, nfit)
            else:
                nc = nfit
            gamma = 2.38 / sqrt(2.0 * nc * nav)
            if uniform() < 0.1: # 10% of the dreamz steps, hop between modes
                gamma = 1.0
            '''
            Make a list of Z index values at the same temperature as X[i].
            '''
            thisTemp = XT[i]
            '''
            Sample indices at random. Take into account averaging.
            '''
            r = sampleIndexWithoutReplacement(nn,2*nav)*nChains + i
            Z1 = sum(Z[r[0:nav],:], axis=0)
            Z2 = sum(Z[r[nav:],:], axis=0)
            '''
            Create an evolution vector from those samples.
            '''
            delta = Z2 - Z1
            '''
            Differential evolution step: make a trial, new solution.
            '''
            x = copy(X[i])
            pp = copy(parIndex)
            if nc  > 0 and nc < ndim: # deal with limited number of crossovers
                pp = sampleWithoutReplacement(parIndex, nc)
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
                    xp = func(x) + pf
                else:
                    xp = func(x, *args)  + pf
                alpha = (xp - XP[i]) / thisTemp
                dice = log(uniform())
                if dice < alpha:
                    X[i] = x
                    XP[i] = xp 
                    if thisTemp == 1: success += [1]
                else:
                    if thisTemp == 1: success += [0]
            else:
                if thisTemp == 1: success += [0]

    return X, XP, success

'''
Swap step: this function allows low temperature and high temperature
solutions to talk to each other. Take a random parameter vector in X
and try swapping it with the next higher temperature parameter vector
in X. The acceptance ratio is based on Earl and Deem (2008). 
'''
def dreamzPTSwap(func, priorFunc, X, XP, XT, args=None):
    TU = unique(XT)
    nChains, ndim = shape(X)
    for i in range(len(TU)-1,0,-1):
        Tlow = TU[i-1]
        Thigh = TU[i]
        betaLow = 1.0 / Tlow
        betaHigh = 1.0 / Thigh

    # pick chains with matching temperatures to swap
        idx = where(XT == Tlow)[0]
        j = sampleWithoutReplacement(idx, 1)[0]
        Plow = XP[j]
        idx = where(XT == Thigh)[0]
        i = sampleWithoutReplacement(idx, 1)[0]
        Phigh = XP[i]

        alpha = (Phigh - Plow)*(betaLow - betaHigh)
        if log(uniform()) < alpha:
            xx = copy(X[j])
            xxp = copy(XP[j])
            X[j] = X[i]
            XP[j] = XP[i]
            X[i] = xx
            XP[i] = xxp
    return X, XP


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
    if all(fitVector == None):
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
    success = []
    #idx = arange(nHist)
    for i in range(nChains):
        '''
        Sample two indices at random.
        '''
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
            X[i] = x
            XP[i] = xp + pf
            success += [1]
        else:
            success += [0]
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
    if all(fitVector == None):
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
    if all(fitVector == None):
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
    if all(fitVector == None):
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)
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

    if all(fitVector == None):
        parIndex = arange(ndim)
    else:
        parIndex = copy(fitVector)

    d = len(parIndex)
    #idz = arange(M)
    success = []
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
            x[parIndex] = x[parIndex] + diff[parIndex]
            pf = priorFunc(x)
            xp = -inf
            if isfinite(pf):
                if args == None:
                        xp = func(x) + pf
                else:
                    xp = func(x, *args) + pf
                bias = pow(vnorm((x-z), 2) / \
                    vnorm((X[i] - z), 2), d-1)
                if log(uniform()) < xp - XP[i] + log(bias):
                    # success!
                    success += [1]
                    XNEW[i] = copy(x)
                    XPNEW[i] = xp
                else:
                    success += [0]
            else:
                success += [0]
        else:
            success += [0]
    return XNEW, XPNEW, success

