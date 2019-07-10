import numpy as np
import matplotlib.pyplot as plt
from testRho import *


def getX(b,xMin,xMax,bMin,bMax):
    return (b-bMin)/(bMax-bMin) * (xMax-xMin) + xMin

def getB(x,xMin,xMax,bMin,bMax):
    return (x-xMin)/(xMax-xMin) * (bMax-bMin) + bMin


def integrate_G_BetasTrapezoid(betas,rho,t):
    integrand = []
    for i,b in enumerate(betas):
        if b > 0.0: integrand.append(rho[i]*np.sin(b*t)/b)
        else:       integrand.append(rho[i])
    return np.trapz(integrand,x=betas)


def integrate_G_XsTrapezoid(betas,rho,t):
    xs = [getX(b,-1,1,betas[0],betas[-1]) for b in betas]
    integrand = []
    for i,x in enumerate(xs):
        b = getB(x,-1,1,betas[0],betas[-1])
        integrationThing = (betas[-1]-betas[0])/(2.0)
        if b > 0.0:
            integrand.append(rho[i]*np.sin(b*t)/b*integrationThing)
        else:
            integrand.append(rho[i]/integrationThing)
    return np.trapz(integrand,x=xs)

def integrate_G_XsGaussLegendre(betas,rho,t,N):
    xs = [getX(b,-1,1,betas[0],betas[-1]) for b in betas]
    integrationVal = 0.0
    points,weights = np.polynomial.legendre.leggauss(N)
    for i in range(N):
        x = points[i]
        w = weights[i]
        b = getB(x,-1,1,betas[0],betas[-1])
        integrationThing = (betas[-1]-betas[0])/(2.0)
        rhoVal = getVal(xs,x,rho)
        if b > 1e-6:
            integrationVal += (w*rhoVal/b*np.sin(b*t)*integrationThing)
        else:
            integrationVal += (w*rhoVal*integrationThing)
    return integrationVal


def getG(betas,rho,t,N):
    return integrate_G_XsGaussLegendre(betas,rho,t,N)

if __name__=='__main__':

    kbT = 0.025
    betas = [x/kbT for x in energies[:]]
    invArea = 1.0/np.trapz(rho,x=betas)
    rho = [invArea*x for x in rho]
    t = 1.0

    print("\nintegrating b1->bmax db'")
    print(integrate_G_BetasTrapezoid(betas,rho,t))

    print("\nintegrating -1->1 dx")
    print(integrate_G_XsTrapezoid(betas,rho,t))

    print("\nintegrating with gauss-legendre")
    print(integrate_G_XsGaussLegendre(betas,rho,t,100))
    print()

    toPrint = False
    if toPrint:
        integrationVals = [integrate_G_XsGaussLegendre(betas,rho,t,N) for N in [1,5,10,20,50,100,500,1000]]
        plt.plot(integrationVals)
        x1 = integrate_G_XsTrapezoid(betas,rho,t)
        plt.plot([x1 for i in range(len(integrationVals))])
        plt.show()

















