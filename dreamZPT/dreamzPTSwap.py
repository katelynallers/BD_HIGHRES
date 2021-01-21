'''
Jack Gallimore & Maggie Overstreet, Bucknell University, 2015
'''
'''
Swap step: this function allows low temperature and high temperature
solutions to talk to each other. Take a random parameter vector in X
and try swapping it with the next higher temperature parameter vector
in X. The acceptance ratio is based on Earl and Deem (2008). 
'''
def dreamzPTSwap(func, priorFunc, X, XP, XT):
    nChains, ndim = shape(X)
    XNEW = copy(X)
    XPNEW = copy(XP)
    i = sampleWithoutReplacement(range(1, nChains), 1)[0]
    j = i-1
    betaLow = 1.0 / XT[j]
    betaHigh = 1.0 / XT[i]
    pfi = priorFunc(X[i])
    pfj = priorFunc(X[j])
    alpha = exp((XP[i] - XP[j])*(betaHigh - betaLow) + pfi - pfj)
    if uniform() < alpha:
        XNEW[i] = X[j]
        XNEW[j] = X[i]
        XPNEW[i] = XP[j]
        XPNEW[j] = XP[i]
    return XNEW, XPNEW

