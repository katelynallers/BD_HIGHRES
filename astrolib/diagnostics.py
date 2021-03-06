'''
Jack Gallimore, Bucknell University, 2015
'''
from numpy import *

# from http://stackoverflow.com/questions/16044491/statistical-scaling-of-autocorrelation-using-numpy-fft
def autocorrelation(x):
    """
    Compute autocorrelation using FFT
    """
    x = asarray(x)
    N = len(x)
    x = x-x.mean()
    s = fft.fft(x, N*2-1)
    result = real(fft.ifft(s * conjugate(s), N*2-1))
    result = result[:N]
    result /= result[0]
    return result

def normalizedACF(p, length=10): # length = window for acf calculation
    xacf = zeros(length-1)
    mp = mean(p)
    mp2 = mean(p*p)
    norm = mp2 - mp*mp
    for i in range(1,length):
        xacf[i-1] = (mean(p[i:]*p[:-i]) - mp*mp) / norm
    return xacf
    

def fullACF(p, length=100): # length = window for acf calculation
    xacf = ones(length)
    mp = mean(p)
    mp2 = mean(p*p)
    norm = mp2 - mp*mp
    for i in range(1,length):
        xacf[i] = (mean(p[i:]*p[:-i]) - mp*mp) / norm
    return xacf


    

# Estimate integrated autocorrelation time.
# Sokal's (1989) technique: calculate truncated sum (integration) 
# within a window. Increase window length T 
# until T >= c * IAT, c = cutoff level. 
# Smaller cutoff introduces bias in estimate, but filters noise.
# Larger cutoff reduces bias in estimate, but increases noise.
# 
# Note: remove burn-in before applying this function.
#
"""
old version
def truncatedIAT(p, cutoff=4): # p = single parameter vector
    crit = inf
    T = 1
    while(T < crit):
        T += 1
        # Most authors seem to prefer this form for iat.
        iat = 1.0 + 2.0*sum(normalizedACF(p, length=T))
        # Note: to get Sokal's IAT, multiply by 0.5. 
        crit = cutoff * iat
    return iat
"""
def truncatedIAT(p, cutoff=4): # p = single parameter vector
    crit = inf
    T = 0
    #xacf = fullACF(p)
    xacf = autocorrelation(p)[1:] # this is a little faster(?) and cleaner
    N = len(p)
    iat = 1.0
    while(T < crit and T < N):
        T += 1
        # Most authors seem to prefer this form for iat.
        iat += 2.0*xacf[T-1]
        # Note: to get Sokal's IAT, multiply by 0.5. 
        crit = cutoff * iat
    return iat

from scipy.optimize import curve_fit
decay = lambda x, tau: exp(-x/tau)
def fittedIAT(p): # p = single parameter vector
    # assumes an exponential decay form for autocorrelation
    crit = inf
    T = 0
    # first, determine the iterative decay iat. There should be a faster
    # way to do this.
    iat, acf = iterativeDecayIATwithACF(p)
    print 'iterativeDecayIAT: ', iat
    #xacf = fullACF(p, length=2*temp) # this is a little faster(?) and cleaner
    lag = arange(len(xacf))
    tau, unc = curve_fit(decay, lag, xacf)
    # Stabilize things by interpolation?
    print 'fitted IAT: ', 1.0 + 2.0 * tau[0]
    
    return 1.0 + 2.0 *tau[0]

def singleLag(p, lag):
    mp = mean(p)
    mp2 = mean(p*p)
    norm = mp2 - mp*mp
    if lag > 0:
        xacf = (mean(p[lag:]*p[:-lag]) - mp*mp) / norm
    else:
        xacf = 1.0
    return xacf

def iterativeDecayIATwithACF(p, lagMax=100000, threshold=0.368):
    acf = []
    n = len(p)
    if n > lagMax: # a little sanity
        p = p[-lagMax:]
        n = lagMax
    lag = 1
    val = 1.0
    acf += [val]
    while lag < n and val > threshold: # look for exp(-1.0)
        oldVal = val
        val = singleLag(p, lag)
        lag += 1
        acf += [val]
    lag = float(lag - 1) + (threshold - oldVal) / (val - oldVal)
    iat = 1.0 + 2.0 * lag
    return iat, acf

def iterativeDecayIAT(p, lagMax=100000, threshold=0.368):
    n = len(p)
    if n > lagMax: # a little sanity
        p = p[-lagMax:]
        n = lagMax
    lag = 1
    val = 1.0
    while lag < n and val > threshold: # look for exp(-1.0)
        oldVal = val
        val = singleLag(p, lag)
        lag += 1
    lag = float(lag - 1) + (threshold - oldVal) / (val - oldVal)
    iat = 1.0 + 2.0 * lag
    return iat
    

def avgSquaredJumpDistance(Z): # Z = history array, filtered by temperature
    # Best to filter out burn-in before calling this function!
    dZ2 = (Z[1:] - Z[:-1])**2 # step sizes squared
    d2 = sum(dZ2, axis=1) # euclidean step distances, squared
    return mean(d2)


def autoburn(chisq, nintervals):
    '''
    Given a chisq trace, automatically determine burn-in with a
    geweke-like statistic.
    '''
    n = len(chisq)
    nn = n / nintervals
    nburnfinal = None
    ### calc Geweke-like score - use end as reference
    mm = mean(chisq[-nn:])
    ss = std(chisq[-nn:])
    # iteratively refine autoburn interval
    while(nn > 20): # at least 20 samples for statistics
        g = inf
        i = 0
        while(i < nintervals-1 and g > 1.2):
            nburn = i*nn
            n1 = (i+1)*nn
            m0 = mean(chisq[nburn:n1])
            g = abs(m0 - mm) / ss
            i += 1
        if g < 1.2:
            nburnfinal = nburn
        n = nburn
        nn = n / nintervals
    if nburnfinal == None:
        print 'Not converged'
        nburnfinal = len(chisq)/2
        print 'Returning burn = len/2 = ', nburnfinal
    return nburnfinal

        
        
    
