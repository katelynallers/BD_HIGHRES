'''
Wrappers for the DREAMZ stepper. Meant to be a nearly drop-in
replacement for other optimizers.
'''
from stepper import dreamzStep, dreamzPTStep, snookerStep, dreamzPTSwap
import numpy as np
import time
import math

def round_sigfigs(num, sig_figs):
    """Round to specified number of sigfigs.

    >>> round_sigfigs(0, sig_figs=4)
    0
    >>> int(round_sigfigs(12345, sig_figs=2))
    12000
    >>> int(round_sigfigs(-12345, sig_figs=2))
    -12000
    >>> int(round_sigfigs(1, sig_figs=2))
    1
    >>> '{0:.3}'.format(round_sigfigs(3.1415, sig_figs=2))
    '3.1'
    >>> '{0:.3}'.format(round_sigfigs(-3.1415, sig_figs=2))
    '-3.1'
    >>> '{0:.5}'.format(round_sigfigs(0.00098765, sig_figs=2))
    '0.00099'
    >>> '{0:.6}'.format(round_sigfigs(0.00098765, sig_figs=3))
    '0.000988'
    """
    if num != 0:
        return np.round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0


def reportElapsedTime(startTime):
    elapsed = time.time() - startTime
    rElapsed = elapsed
    ts = 'seconds'
    if elapsed > 60.:
        rElapsed = elapsed / 60.
        ts = 'minutes'
    if elapsed > 3600:
        rElapsed = elapsed / 3600.
        ts = 'hours'
    if elapsed > 3600*24.:
        rElapsed = elapsed / 3600. / 24.
        ts = 'days'
    if elapsed > 3600.*24.*365.:
        rElapsed = elapsed / (3600.*24*365.)
        ts = 'years'
    print 'Elapsed time: ', round_sigfigs(rElapsed, 3), ts
    return

def reportTimeRemaining(startTime, i, niter):
    elapsed = time.time() - startTime # in seconds
    j = i
    if j == 0: j = 1
    tpi = elapsed / float(j)
    timeRemaining = tpi * (niter-i)
    tr = timeRemaining
    ts = 'seconds'
    if timeRemaining > 60.:
        tr = timeRemaining / 60.
        ts = 'minutes'
    if timeRemaining > 3600:
        tr = timeRemaining / 3600.
        ts = 'hours'
    if timeRemaining > 3600*24.:
        tr = timeRemaining / 3600. / 24.
        ts = 'days'
    if timeRemaining > 3600.*24.*365.:
        tr = timeRemaining / (3600.*24*365.)
        ts = 'years'
    print 'Time per iteration: ', round_sigfigs(tpi, 3), ' seconds'
    print 'Time Remaining: ', round_sigfigs(tr, 3), ts
    return

def ncrDraw(ncrScl, ncrMax):
    selectFrom = np.arange(ncrMax) + 1
    arg = (1.0*selectFrom) / (1.0 * ncrScl)
    probs = np.exp(-arg*arg)
    cdf = np.cumsum(probs) / np.sum(probs)
    dice = np.random.uniform()
    ncr = selectFrom[np.where(dice < cdf)[0][0]]
    return ncr

def DREAM(X, logl, logp, niter, thin=1, ncr=0, navg=1, fitVector=None, \
              args=None):
    '''
    X = initial guess at parameters; shape = (nchains, ndim)

    logl = function returns log likelihood (usually -chisq / 2)

    logp = function returns log prior probability (usually -inf or 0 to 
    define boundaries).
    '''
    nchains, ndim = np.shape(X)
    XP = []
    for x in X:
        xp = logl(x)  + logp(x)
        XP = XP + [xp]
    XP = np.array(XP)

    # allocate memory for results
    Z = np.zeros((niter*nchains, ndim))
    ZP = np.zeros(niter*nchains)
    Z[0:nchains] = X
    ZP[0:nchains] = XP

    # store record of successes/failures
    successVec = []
    for i in range(niter):
        if i % 1000 == 0: print i, '/', niter
        for j in range(thin):
            nav = np.random.randint(1,navg+1)
            nc = ncr
            gammaZ = 2.38 / np.sqrt(2.0 * nav * nc) # scaling of difference vectors
            if np.random.uniform() < 0.1: # 10% of the time, hop between modes
                gammaZ = 1.0
                nc = ndim
            X, XP, success = dreamzStep(logl, logp, X, XP, X, XP, \
                                            ncr=nc, navg=nav, \
                                            gamma=gammaZ, \
                                            fitVector=fitVector, \
                                            args=args)
            successVec += success
        # update solutions
        i0 = i*nchains
        i1 = i0 + nchains
        Z[i0:i1] = X
        ZP[i0:i1] = XP

    return Z, ZP, successVec


def DREAMZ(ZZ, logl, logp, nchains, niter, thin=1, ncr=0, navg=1, fitVector=None, \
              args=None):
    '''
    ZZ = initial guess at parameters; shape = (ninit, ndim)

    logl = function returns log likelihood (usually -chisq / 2)

    logp = function returns log prior probability (usually -inf or 0 to 
    define boundaries).
    '''
    ninit, ndim = np.shape(ZZ)
    ZZP = []
    for z in ZZ:
        zp = logl(z)  + logp(z)
        ZZP = ZZP + [zp]
    ZZP = np.array(ZZP)

    # allocate memory for results
    Z = np.zeros((ninit+niter*nchains, ndim))
    ZP = np.zeros(ninit+niter*nchains)
    Z[0:ninit] = ZZ
    ZP[0:ninit] = ZZP
    X = np.copy(ZZ[-nchains:])
    XP = np.copy(ZZP[-nchains:])

    # store record of successes/failures
    successVec = []
    nsofar = ninit
    for i in range(niter):
        if i % 1000 == 0: print i, '/', niter
        for j in range(thin):
            nav = np.random.randint(1,navg+1)
            nc = ncr
            if ncr <= 0:
                gammaZ = 2.38 / np.sqrt(2.0 * ndim) # scaling of difference vectors
            else:
                gammaZ = 2.38 / np.sqrt(2.0 * ncr * nav) # scaling of difference vectors
            if np.random.uniform() < 0.1: # 10% of the time, hop between modes
                gammaZ = 1.0
                nc = ndim
            X, XP, success = dreamzStep(logl, logp, X, XP, Z[0:nsofar], ZP[0:nsofar], \
                                            ncr=nc, navg=nav, \
                                            gamma=gammaZ, \
                                            fitVector=fitVector, \
                                            args=args)
            successVec += success
        # update solutions
        i0 = i*nchains
        i1 = i0 + nchains
        Z[i0:i1] = X
        ZP[i0:i1] = XP
        nsofar += nchains

    return Z, ZP, successVec

'''
Wrapper for DREAM(ZS) algorithm.
Added reflections at boundaries if plow and phigh are defined.
Plow = vector of lower boundary values, one for each parameter.
Phigh = ditto for upper boundary values.
'''

def DREAMZS(ZZ, logl, logp, nchains, niter, thin=1,\
            ncr=0, navg=1, \
            fitVector=None, \
            args=None, verbose=True, snookerProb=0.1, \
            plow=None, phigh=None):
    '''
    ZZ = initial guess at parameters; shape = (ninit, ndim)

    logl = function returns log likelihood (usually -chisq / 2)

    logp = function returns log prior probability (usually -inf or 0 to 
    define boundaries).
    '''
    ninit, ndim = np.shape(ZZ)
    ZZP = []
    for z in ZZ:
        zp = logl(z)  + logp(z)
        ZZP = ZZP + [zp]
    ZZP = np.array(ZZP)

    if fitVector != None:
        nfit = len(fitVector)
    else:
        nfit = ndim

    # allocate memory for results
    Z = np.zeros((ninit+niter*nchains, ndim))
    ZP = np.zeros(ninit+niter*nchains)
    Z[0:ninit] = ZZ
    ZP[0:ninit] = ZZP
    X = np.copy(ZZ[-nchains:])
    XP = np.copy(ZZP[-nchains:])

    # store record of successes/failures
    successVec = []
    nsofar = ninit
    startTime = time.time()
    for i in range(niter):
        if verbose:
            if i % 1000 == 0: 
                print ''
                print i, '/', niter 
                print 'XP = ', XP
                sr = sum(successVec[-1000:]) / 1000.
                if i > 0:
                    reportElapsedTime(startTime)
                    reportTimeRemaining(startTime, i, niter)
                    #print 'nsofar: ', nsofar
                    #print ZP[0:nsofar]
                    #chisq = -2.0 * amax(ZP[0:nsofar])
                    #print 'Chisq = ', chisq
                    print 'Chisq = ', round_sigfigs(-2.0*np.amax(ZP[0:nsofar]),4)
                    print 'Success = ', round_sigfigs(sr, 3)
        for j in range(thin):
            if np.random.uniform() < snookerProb: 
                X, XP, success = snookerStep(logl, logp, X, XP, Z[0:nsofar], \
                                                 fitVector=fitVector, \
                                                 args=args)
            else:
                nav = np.random.randint(1,navg+1)
                if ncr > 0: 
                    #nc = np.random.randint(1,ncr+1)
                    nc = ncrDraw(ncr, nfit)
                else:
                    nc = nfit
                gammaZ = 2.38 / np.sqrt(2.0 * nc * nav)
                if np.random.uniform() < 0.1: # 10% of the dreamz steps, hop between modes
                    gammaZ = 1.0
                    nc = ndim
                    nav = 1
                X, XP, success = dreamzStep(logl, logp, X, XP, Z[0:nsofar],\
                                                ZP[0:nsofar], \
                                                ncr=nc, navg=nav, \
                                                gamma=gammaZ, \
                                                fitVector=fitVector, \
                                                args=args, plow=plow, phigh=phigh)
            successVec += success
        # update solutions
        i0 = ninit + i*nchains
        i1 = i0 + nchains
        Z[i0:i1] = np.copy(X)
        ZP[i0:i1] = np.copy(XP)
        #nsofar += nchains
        nsofar = i1

    return Z, ZP, successVec

def DREAMZSA(ZZ, logl, logp, nchains, niter, thin=1, ncr=0, navg=1, \
                fitVector=None, \
                args=None, verbose=True, snookerProb=0.1):
    '''
    ZZ = initial guess at parameters; shape = (ninit, ndim)

    logl = function returns log likelihood (usually -chisq / 2)

    logp = function returns log prior probability (usually -inf or 0 to 
    define boundaries).

    This is an adaptive-gamma version of dreamZS.
    '''
    ninit, ndim = np.shape(ZZ)
    ZZP = []
    for z in ZZ:
        zp = logl(z)  + logp(z)
        ZZP = ZZP + [zp]
    ZZP = np.array(ZZP)

    if fitVector != None:
        nfit = len(fitVector)
    else:
        nfit = ndim

    # allocate memory for results
    Z = np.zeros((ninit+niter*nchains, ndim))
    ZP = np.zeros(ninit+niter*nchains)
    Z[0:ninit] = ZZ
    ZP[0:ninit] = ZZP
    X = np.copy(ZZ[-nchains:])
    XP = np.copy(ZZP[-nchains:])

    # store record of successes/failures
    successVec = []
    nsofar = ninit
    currentGamma = 2.38  # for some reason dividing by two works better
    startTime = time.time()
    for i in range(niter):
        if verbose:
            if i % 1000 == 0: 
                print ''
                print i, '/', niter
                if i > 0:
                    reportElapsedTime(startTime)
                    reportTimeRemaining(startTime, i, niter)
                print 'Chisq = ', round_sigfigs(-2.0*np.amax(ZP[0:nsofar]),4)
                sr = sum(successVec[-1000:]) / 1000.
                # adapt currentGamma to optimize success rate
                cr = 1.0
                if sr < 0.05: 
                    cr = 0.5
                if sr < 0.001:
                    cr = 0.1
                if sr > 0.5:
                    cr = 1.1
                if sr > 0.75:
                    cr = 2.0
                if sr > 0.95:
                    cr = 10.0
                currentGamma *= cr
                if i > 0: print 'Success = ', round_sigfigs(sr, 3)
                print 'Gamma = ', round_sigfigs(currentGamma, 3)
                    
        for j in range(thin):
            if np.random.uniform() < snookerProb: 
                X, XP, success = snookerStep(logl, logp, X, XP, Z[0:nsofar], \
                                                 fitVector=fitVector, \
                                                 args=args)
            else:
                nav = np.random.randint(1,navg+1)
                if ncr > 0: 
                    #nc = np.random.randint(1,ncr+1)
                    nc = ncrDraw(ncr, nfit)
                else:
                    nc = nfit
                gammaZ = currentGamma / np.sqrt(2.0 * nc * nav) 
                if np.random.uniform() < 0.1: # 10% of the dreamz steps, hop between modes
                    gammaZ = 1.0
                    nc = ndim
                    nav = 1
                X, XP, success = dreamzStep(logl, logp, X, XP, Z[0:nsofar],\
                                                ZP[0:nsofar], \
                                                ncr=nc, navg=nav, \
                                                gamma=gammaZ, \
                                                fitVector=fitVector, \
                                                args=args)
            successVec += success
        # update solutions
        i0 = ninit + i*nchains
        i1 = i0 + nchains
        Z[i0:i1] = X
        ZP[i0:i1] = XP
        #nsofar += nchains
        nsofar = i1

    return Z, ZP, successVec


def DREAMZPT(ZZ, logl, logp, XT, niter, thin=1, ncr=0, navg=1, \
                fitVector=None, \
                args=None, verbose=True,plow=None, phigh=None):
    '''
    ZZ = initial guess at parameters; shape = (ninit, ndim)

    logl = function returns log likelihood (usually -chisq / 2)

    logp = function returns log prior probability (usually -inf or 0 to 
    define boundaries).
    '''
    ninit, ndim = np.shape(ZZ)
    ZZP = []
    for z in ZZ:
        zp = logl(z)  + logp(z)
        ZZP = ZZP + [zp]
    ZZP = np.array(ZZP)

    nchains = len(XT)
    ZT = np.zeros(niter*nchains+ninit) # temperature array
    
    # allocate memory for results
    Z = np.zeros((ninit+niter*nchains, ndim))
    ZP = np.zeros(ninit+niter*nchains)
    ZT = np.zeros(ninit+niter*nchains)
    Z[0:ninit] = ZZ
    ZP[0:ninit] = ZZP
    X = np.copy(ZZ[-nchains:])
    XP = np.copy(ZZP[-nchains:])
    ZT[0:ninit] = np.tile(XT, ninit / nchains + 1)[0:ninit]

    # store record of successes/failures
    successVec = []
    nsofar = ninit
    startTime = time.time()
    for i in range(niter):
        if verbose:
            if i % 1000 == 0 and i != 0: 
                print ''
                print i, '/', niter
                print 'Chisq = ', round_sigfigs(-2.0*np.amax(ZP[0:nsofar]),4)
                if len(successVec) > 0:
                    sr = sum(successVec) / (1.0 * len(successVec))
                    print 'Success = ', round_sigfigs(sr, 3)
                reportElapsedTime(startTime)
                reportTimeRemaining(startTime, i, niter)
                successVec = []
        X, XP, success = dreamzPTStep(logl, logp, X, XP, XT, Z[0:nsofar],\
                                            ZP[0:nsofar], ZT[0:nsofar], \
                                            ncr=ncr, navg=navg, \
                                            fitVector=fitVector, \
                                            args=args, plow=None, phigh=None, \
                                            thin=thin)
        successVec += success
        X, XP = dreamzPTSwap(logl, logp, X, XP, XT, args=args)
        # update solutions
        i0 = ninit + i*nchains
        i1 = i0 + nchains
        Z[i0:i1] = np.copy(X)
        ZP[i0:i1] = np.copy(XP)
        ZT[i0:i1] = np.copy(XT)
        #nsofar += nchains
        nsofar = i1

    return Z, ZP, ZT, successVec

