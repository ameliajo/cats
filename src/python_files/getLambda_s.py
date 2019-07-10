import numpy as np
import matplotlib.pyplot as plt
from testRho import *

def getX(b,xMin,xMax,bMin,bMax):
    return (b-bMin)/(bMax-bMin) * (xMax-xMin) + xMin

def getB(x,xMin,xMax,bMin,bMax):
    return (x-xMin)/(xMax-xMin) * (bMax-bMin) + bMin



def integrate_Lambda_BetasTrapezoid(betas,rho):
    integrand = []
    for i,b in enumerate(betas):
        integrand.append(rho[i]/(b*np.tanh(b*0.5)))
    return np.trapz(integrand,x=betas)

def integrate_Lambda_XsTrapezoid(betas,rho):
    xs = [getX(b,-1,1,betas[0],betas[-1]) for b in betas]
    integrand = []
    for i,x in enumerate(xs):
        b = getB(x,-1,1,betas[0],betas[-1])
        integrationThing = (betas[-1]-betas[0])/(2.0)
        integrand.append(rho[i]/(b*np.tanh(b*0.5))*integrationThing)
    return (np.trapz(integrand,x=xs))

def integrate_Lambda_XsGaussLegendre(betas,rho,N):
    xs = [getX(b,-1,1,betas[0],betas[-1]) for b in betas]
    integrationVal = 0.0
    points,weights = np.polynomial.legendre.leggauss(N)
    for i in range(N):
        x = points[i]
        w = weights[i]
        b = getB(x,-1,1,betas[0],betas[-1])
        integrationThing = (betas[-1]-betas[0])/(2.0)
        rhoVal = getVal(xs,x,rho)
        integrationVal += (w*rhoVal/(b*np.tanh(b*0.5))*integrationThing)
    return integrationVal


def getSmallContrib(betas,rho):
    c = rho[1]/(betas[1]*betas[1])
    return 2.0*c*betas[1] + c*betas[1]**3/18.0 



if __name__=='__main__':

    kbT = 8.6173332e-5*296.0
    betas = [x/kbT for x in energies[:]]
    invArea = 1.0/np.trapz(rho,x=betas)
    print("area",1.0/invArea)
    rho = [invArea*x for x in rho]
    print(rho[0],rho[1])

    print("\nintegrating b1->bmax db'")
    print(integrate_Lambda_BetasTrapezoid(betas[1:],rho))
    
    print("\nintegrating -1->1 dx")
    print(integrate_Lambda_XsTrapezoid(betas[1:],rho))
    x1 = integrate_Lambda_XsTrapezoid(betas[1:],rho)
    
    print("\nintegrating with gauss-legendre")
    print(integrate_Lambda_XsGaussLegendre(betas[1:],rho,100))
    toPrint = False
    if toPrint:
        integrationVals = [integrate_Lambda_XsGaussLegendre(betas[1:],rho,N) for N in [1,5,10,20,50,100,500,1000]]
        plt.plot(integrationVals)
        x1 = integrate_Lambda_XsTrapezoid(betas[1:],rho)
        plt.plot([x1 for i in range(len(integrationVals))])
        plt.show()

    print("\nsmall contribution")
    print(getSmallContrib(betas,rho))























