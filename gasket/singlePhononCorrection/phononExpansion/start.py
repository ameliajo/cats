from getDebyeWaller import getDebyeWaller
from normalize      import normalize 
from math import sinh, exp



def start(betas, rho, wgt):
    P = [rho[1]/betas[1]**2] + \
        [rho[b]/(2.0*betas[b]*sinh(betas[b]*0.5)) for b in range(1,len(rho))]
    P = normalize(betas,P,wgt)
    lambda_s = getDebyeWaller(betas,P)
    T1 = [P[b]*exp(betas[b]*0.5)/lambda_s for b in range(len(betas))]
    return lambda_s, T1
