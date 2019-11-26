from math import sinh
import numpy as np


def normalize(betas,P,wgt):
    integrand = [P[b]*2.0*betas[b]*sinh(betas[b]*0.5) for b in range(len(P))]
    invSum = wgt/np.trapz(integrand,betas)
    P = [Pval * invSum for Pval in P]
    return P 



