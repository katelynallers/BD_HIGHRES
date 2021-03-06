'''
Jack Gallimore & Maggie Overstreet, Bucknell University, 2015
'''
'''
Perform differential evolution (DREAMZ). Parallel tempering
modification: only allow single temperature evolution. In other words,
a given X[i] only is modified by Z values at the same temperature.
'''
def dreamzPTStep(func, priorFunc, \
                     X, XP, XT, Z, ZP, ZT, \
                     gamma=1.0, eps=1.e-6, e=0.05, \
                     fixed=None):
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
        if fixed != None:
            x[fixed] = X[fixed]
        xp = func(x)
        '''
        Acceptance ratio setp. Notice the inclusion of temperature.
        '''
        pfi = priorFunc(X[i])
        pf = priorFunc(x)
        alpha = exp((xp - XP[i])/XT[i] \
                        + pf - pfi)
        if uniform() < alpha:
            X[i] = x
            XP[i] = xp
    return X, XP

