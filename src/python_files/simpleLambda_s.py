import numpy as np
from testRho import *
from phononDists import *

def getVal(energies,val,vec):
    for i in range(len(energies)-1):
        if energies[i] <= val and val <= energies[i+1]:
            m = (vec[i+1]-vec[i])/(energies[i+1]-energies[i])
            b = vec[i]
            return m*(val-energies[i])+b
    return 0.0
 
doQuad = True
doQuad = False

if doQuad:
    integrationVal = 0.0
    N = 100
    points,weights = np.polynomial.legendre.leggauss(N)
    xs = [2*b/betas[-1]-1 for b in betas]
    for i in range(N):
        x, w = points[i], weights[i]
        b = betas[-1]*(x+1)*0.5
        integrationThing = betas[-1]*0.5
        rhoVal = getVal(xs,x,rho)
        integrationVal += w*rhoVal*integrationThing
    rho = [x/integrationVal for x in rho]
    print("\nintegration val",integrationVal,"\n")
else:
    invArea = 1.0/np.trapz(rho,x=betas)
    rho = [invArea*x for x in rho]

################################################################################
################################################################################

P = [rho[1]/(betas[1]*betas[1])] + [rho[i]/(2.0*betas[i]*np.sinh(betas[i]*0.5)) for i in range(1,len(rho))]
fullP = P[::-1] + P[1:]
fullBetas = [-b for b in betas][::-1] + betas[1:]

################################################################################
################################################################################

if doQuad:
    xs = [b/betas[-1] for b in fullBetas]
    integrationVal = 0.0
    N = 10
    points,weights = np.polynomial.legendre.leggauss(N)
    for i in range(N):
        x, w = points[i], weights[i]
        b = betas[-1]*x
        integrationThing = betas[-1]
        pVal = getVal(xs,x,fullP)
        print("------- ",pVal*np.exp(-b*0.5)*integrationThing)
        integrationVal += w*pVal*np.exp(-b*0.5)*integrationThing
    print(integrationVal)
else:
    integrand = [fullP[i]*np.exp(-fullBetas[i]*0.5) for i in range(len(fullP))]
    print(np.trapz(integrand,x=fullBetas))







