from math import cosh
import numpy as np


def getDebyeWaller(betas,P):

    integrand = [P[b]*2.0*cosh(betas[b]*0.5) for b in range(len(betas))]
    return np.trapz(integrand,betas)



    
