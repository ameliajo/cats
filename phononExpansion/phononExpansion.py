import numpy as np
from numpy import sinh, cosh, exp
from convolution import *


def interpolate(vec, x, delta):
    if x < 0 or x > (len(vec)-1)*delta:
        return 0.0
    for i in range(len(vec)-1):
        if i*delta <= x and (i+1)*delta > x:
            b = vec[i]
            m = (vec[i+1]-vec[i])/delta
            return m*(x-i*delta)+b
    return 0.0

def contin(nphon, delta, rho, alpha, beta):
    betaGrid = [delta*i for i in range(len(rho))]
    P = [rho[1]/betaGrid[1]**2] + \
        [rho[i]/(2*betaGrid[i]*sinh(betaGrid[i]*0.5)) for i in range(1,len(rho))]
    integrand = [P[i]*2.0*sinh(betaGrid[i]*0.5) for i in range(len(rho))]
    invSum = 1.0/np.trapz(integrand,x=betaGrid)
    P = [invSum*val for val in P]
    lambda_s = np.trapz([2.0*P[i]*cosh(betaGrid[i]*0.5) for i in range(len(rho))],x=betaGrid)
    t1 = [P[b]*exp(betaGrid[b]*0.5)/lambda_s for b in range(len(P))]

    sab = [0.0]*len(alpha)*len(beta)
    xa = [1.0]*len(alpha)

    tnow  = [0.0]*nphon*len(t1)
    tlast = [0.0]*nphon*len(t1)
    for i in range(len(t1)):
        tnow[i]  = t1[i]
        tlast[i] = t1[i]

    nNext, nLast = len(t1), len(t1)

    lambda_alpha = [lambda_s*alphaVal for alphaVal in alpha]
    exp_lambda_alpha = [exp(-lambda_alpha_val) for lambda_alpha_val in lambda_alpha]

    for n in range(nphon):
        if (n>0):
            nNext = len(t1)+nLast-1
            tnow = convol(t1,tlast,delta,nNext)
        inv_n = 1.0/(n+1)
        for a in range(len(alpha)):
            xa[a] *= lambda_alpha[a]*inv_n
            exx    = exp_lambda_alpha[a]*xa[a]

            for b in range(len(beta)):
                add = exx * interpolate(tnow, beta[b], delta)
                if add > 1e-30:
                    sab[b+a*len(beta)] += add
        if n > 0:
            for i in range(nNext):
                tlast[i] = tnow[i]
            nLast = nNext
    return sab



if __name__=='__main__':
    from dos import *
    import matplotlib.pyplot as plt
    print("This probably has to be debugged, doesnt match NJOY")

    alpha = np.linspace(0,5,20)
    beta  = list(np.linspace(0,1,30)) + list(np.linspace(1.1,5,30)) + list(np.linspace(5.1,20,20))

    delta = delta_water
    rho = rho_water[:]
    betaGrid = [delta*i for i in range(len(rho))]
    sab = contin(20,delta,rho,alpha,beta)

    for a in range(len(alpha)):
        sabChunk = [sab[b+a*len(beta)] for b in range(len(beta))]
        plt.plot(beta,sabChunk)
    plt.show()


