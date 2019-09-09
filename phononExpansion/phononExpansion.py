import numpy as np
from numpy import sinh, exp, cosh


def getConvolAtPoint(i, delta, t1, t2):
    sumVal = 0.0
    for j in range(-len(t1)+1,len(t1)):
        if ( i - j >= len(t2) ):
            return sumVal
        expVal = exp(j*delta) if j < 0 else exp((i-j)*delta) if i < j else 1.0
        toAdd = t1[abs(j)]*t2[abs(i-j)]*expVal
        sumVal += 0.5*toAdd if j == -len(t1)+1 or j == len(t1)-1 else toAdd
    return sumVal


def convol(t1, t2, delta, nn):
    t3 = [0.0]*nn
    for i in range(0, nn):
        t3[i] = getConvolAtPoint(i,delta,t1,t2) * delta
        t3[i] = 0 if t3[i] < 1e-30 else t3[i]
    return t3








def contin(nphon, delta, rho, alpha, beta):
    betaGrid = [delta*i for i in range(len(rho))]
    P = [rho[1]/betaGrid[1]**2] + \
        [rho[i]/(2*betaGrid[i]*sinh(betaGrid[i]*0.5)) for i in range(1,len(rho))]
    integrand = [P[i]*2*betaGrid[i]*sinh(betaGrid[i]*0.5) for i in range(len(P))]
    invArea = 1.0/np.trapz(integrand,x=betaGrid)
    P = [invArea * val for val in P]
    lambda_s = np.trapz([P[i]*2*cosh(-betaGrid[i]*0.5) for i in range(len(P))],\
                     x=betaGrid)
    # Remember this is actually Tn(-b) not Tn(b), so this is multiplied by a
    # positive 0.5 exponential rather than a negative.
    T1 = [P[b]*exp(betaGrid[b]*0.5)/lambda_s for b in range(len(P))]

    nNext = len(t1)
    nLast = len(t1)
    Tlast = [0.0]*nphon*len(t1)
    TNext = [0.0]*nphon*len(t1)
    for i in range(len(T1)):
        Tlast[i] = T1[i]
        TNext[i] = T1[i]

    for n in range(nphon):
        if n > 0:
            nNext = len(t1)+len(tlast)-1
            tnow = convol(T1,Tlast,delta,nNext)
        inv_n = 1.0/(n+1)
        for a in range(len(alpha)):


    return lambda_s





















