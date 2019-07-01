import numpy as np
import matplotlib.pyplot as plt
from testRho import *

def getX(b,xMin,xMax,bMin,bMax):
    return (b-bMin)/(bMax-bMin) * (xMax-xMin) + xMin

def getB(x,xMin,xMax,bMin,bMax):
    return (x-xMin)/(xMax-xMin) * (bMax-bMin) + bMin


def integrate_F_BetasTrapezoid(betas,rho,t):
    integrand = []
    for i,b in enumerate(betas):
        integrand.append(rho[i]/(b*np.tanh(b*0.5))*np.cos(b*t))
    return np.trapz(integrand,x=betas)

def integrate_F_XsTrapezoid(betas,rho,t):
    xs = [getX(b,-1,1,betas[0],betas[-1]) for b in betas]
    integrand = []
    for i,x in enumerate(xs):
        b = getB(x,-1,1,betas[0],betas[-1])
        integrationThing = (betas[-1]-betas[0])/(2.0)
        integrand.append(rho[i]/(b*np.tanh(b*0.5))*np.cos(b*t)*integrationThing)
    return (np.trapz(integrand,x=xs))

def integrate_F_XsGaussLegendre(betas,rho,t,N):
    xs = [getX(b,-1,1,betas[0],betas[-1]) for b in betas]
    integrationVal = 0.0
    points,weights = np.polynomial.legendre.leggauss(N)
    for i in range(N):
        x = points[i]
        w = weights[i]
        b = getB(x,-1,1,betas[0],betas[-1])
        integrationThing = (betas[-1]-betas[0])/(2.0)
        rhoVal = getVal(xs,x,rho)
        integrationVal += (w*rhoVal/(b*np.tanh(b*0.5))*np.cos(b*t)*integrationThing)
    return integrationVal


def getSmallContrib(beta1,rho1,t):
    c = rho1/(beta1*beta1)
    return 2.0*c*beta1 + (c/18.0 + c*t*t/3.0)*beta1**3 - (c*t*t/60.0)*beta1**5


def getF(betas,rho,t,N):
    return integrate_F_XsGaussLegendre(betas,rho,t,N) 

if __name__=='__main__':

    kbT = 0.025
    betas = [x/kbT for x in energies[:]]
    invArea = 1.0/np.trapz(rho,x=betas)
    rho = [invArea*x for x in rho]

    t = 1.0

    print("\nintegrating b1->bmax db'")
    print(integrate_F_BetasTrapezoid(betas[1:],rho,t))
    
    print("\nintegrating -1->1 dx")
    print(integrate_F_XsTrapezoid(betas[1:],rho,t))
    x1 = integrate_F_XsTrapezoid(betas[1:],rho,t)
    
    print("\nintegrating with gauss-legendre")
    print(integrate_F_XsGaussLegendre(betas[1:],rho,t,100))
    toPrint = False
    if toPrint:
        integrationVals = [integrate_F_XsGaussLegendre(betas,rho,t,N) for N in [1,5,10,20,50,100,500,1000]]
        plt.plot(integrationVals)
        x1 = integrate_F_XsTrapezoid(betas,rho,t)
        plt.plot([x1 for i in range(len(integrationVals))])
        plt.show()

    print("\nsmall contribution")
    print(getSmallContrib(betas[1],rho[1],t))























