'''
Jack Gallimore, Bucknell University, 2015
'''
'''
Wrappers for the DREAMZ stepper. Meant to be a nearly drop-in
replacement for other optimizers.
'''
from stepper import dreamzStep, snookerStep
from gwstepper import gwStepZ
import numpy as np

def ncrDraw(ncrScl, ncrMax):
    ncrMin = 1
    selectFrom = np.arange(ncrMax) + 1
    arg = (1.0*selectFrom) / (1.0 * ncrScl)
    probs = np.exp(-arg*arg)
    cdf = np.cumsum(probs) / np.sum(probs)
    dice = np.random.uniform()
    ncr = selectFrom[np.where(dice < cdf)[0][0]]
    return ncr
    

def DREAMZHYBRID(ZZ, logl, logp, nchains, niter, thin=1, ncr=0, navg=1, \
                fitVector=None, \
                args=None, scale=2, verbose=True):
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
    X = ZZ[-nchains:]
    XP = ZZP[-nchains:]

    # store record of successes/failures
    successVec = []
    nsofar = ninit
    for i in range(niter):
        if verbose:
            if i % 1000 == 0: 
                print i, '/', niter
                print 'max log P = ', np.amax(ZP[0:nsofar])
                print 'Success = ', np.sum(successVec[-1000:]) / 1000.
        for j in range(thin):
            if np.random.uniform() < 0.5: # try a Goodman Weare Z step
                X, XP, success = gwStepZ(logl, logp, X, XP, Z[0:nsofar],\
                                             fitVector=fitVector, \
                                             scale=scale, \
                                             args=args)
                
            else:
                nav = np.random.randint(1,navg+1)
                if ncr > 0: 
                    #nc = np.random.randint(1,ncr+1)
                    nc = ncrDraw(ncr, nfit)
                else:
                    nc = nfit
                gammaZ = 2.38 / np.sqrt(2.0 * nc * nav) 
                if np.random.uniform() < 0.1: # 10% of the time, hop between modes
                    gammaZ = 1.0
                    '''
                    nc = ndim
                    nav = 1
                    '''
                X, XP, success = dreamzStep(logl, logp, X, XP, Z[0:nsofar],\
                                                ZP[0:nsofar], \
                                                ncr=nc, navg=nav, \
                                                gamma=gammaZ, \
                                                fitVector=fitVector, \
                                                args=args)
                successVec += success
        # every thin-th iteration, try a snooker step
        X, XP, success = snookerStep(logl, logp, X, XP, Z[0:nsofar], \
                                         fitVector=fitVector, \
                                         args=args)
        # update solutions
        i0 = ninit + i*nchains
        i1 = i0 + nchains
        Z[i0:i1] = X
        ZP[i0:i1] = XP
        #nsofar += nchains
        nsofar = i1

    return Z, ZP, successVec

