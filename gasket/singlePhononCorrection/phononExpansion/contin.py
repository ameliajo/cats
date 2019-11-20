from start import *
from convol import *
import numpy as np
from math import exp

def interpolate(yVec, xVal, delta):
    if xVal >= len(yVec)*delta or xVal < 0: return 0.0
    for i in range(len(yVec)-1):
        if i*delta <= xVal and xVal < (i+1)*delta:
            return yVec[i] + (xVal-delta*i) * (yVec[i+1]-yVec[i]) / delta
    return 0.0

def contin(nphon, delta, wgt, rho, alphas, betas):
    betaGrid     = [delta*i for i in range(len(rho))]
    lambda_s, t1 = start(betaGrid,rho,wgt)
    sab = [0.0]*len(alphas)*len(betas)
    xa  = [1.0]*len(alphas)
    tnow  = [0.0]*(len(t1)*nphon)
    tlast = [0.0]*(len(t1)*nphon)
    nLast = len(t1)
    for i in range(len(t1)):
        tnow[i]  = t1[i]
        tlast[i] = t1[i]
    lambda_alpha = [lambda_s*alpha for alpha in alphas]
    exp_lambda_alpha = [exp(-lambda_alpha_val) for lambda_alpha_val in lambda_alpha]

    for n in range(nphon):
        if n > 0:
            nNext = len(t1)+nLast-1
            tnow = convol(t1, tlast, delta, nNext)
        inv_n = 1.0/(n+1)
        for a in range(len(alphas)):
            xa[a] *= lambda_alpha[a] * inv_n
            exx    = exp_lambda_alpha[a]*xa[a]

            for b in range(len(betas)):
                add = exx * interpolate(tnow,betas[b],delta)
                if add > 1e-30:
                    sab[b+a*len(betas)] += add
        if n == 0:
            continue
        for i in range(nNext):
            tlast[i] = tnow[i]
        nLast = nNext
    return sab
