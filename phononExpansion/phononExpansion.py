import numpy as np
import matplotlib.pyplot as plt
from numpy import sinh, cosh, exp


def getLambda_s(rho,betaGrid):
    P = [rho[1]/betaGrid[1]**2] + \
        [rho[i]/(2*betaGrid[i]*sinh(betaGrid[i]*0.5)) for i in range(1,len(rho))]

    integrand = [P[i]*2.0*sinh(betaGrid[i]*0.5) for i in range(len(rho))]
    invSum = 1.0/np.trapz(integrand,x=betaGrid)
    P = [invSum*val for val in P]
    lambda_s = np.trapz([2.0*P[i]*cosh(betaGrid[i]*0.5) for i in range(len(rho))],x=betaGrid)
    return lambda_s




def contin(nphon, delta, rho, alpha, beta):
    betaGrid = [delta*i for i in range(len(rho))]


#delta = 0.03
#rho = [-(i-20)**2+400 for i in range(40)] + [10 for i in range(40,47)] + [-(i-70)**2+500 for i in ra#nge(48,93)] + [0.0]
#betaGrid = [delta*i for i in range(len(rho))]
#plt.plot(rho)
#plt.show()




