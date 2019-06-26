import numpy as np
import matplotlib.pyplot as plt
from testRho import *

bMin=0.1
bMax=10.0

xMin=-1.0
xMax= 1.0

def getX(b,xMin,xMax,bMin,bMax):
    return (b-bMin)/(bMax-bMin) * (xMax-xMin) + xMin

def getB(x,xMin,xMax,bMin,bMax):
    return (x-xMin)/(xMax-xMin) * (bMax-bMin) + bMin


kbT = 0.025

bs = [x/kbT for x in energies[1:]]
bMin, bMax = bs[0], bs[-1]
xs = [getX(b,xMin,xMax,bMin,bMax) for b in bs]

t = 1.0

integrand = []
for i,b in enumerate(bs):
    integrand.append(rho[i]/(b*np.tanh(b*0.5))*np.cos(b*t))
print(np.trapz(integrand,x=bs))


integrand = []
for i,x in enumerate(xs):
    b = getB(x,xMin,xMax,bMin,bMax)
    integrationThing = (bMax-bMin)/(xMax-xMin)
    integrand.append(rho[i]/(b*np.tanh(b*0.5))*np.cos(b*t)*integrationThing)
print(np.trapz(integrand,x=xs))
x1 = (np.trapz(integrand,x=xs))


integrationVals = []
for N in [1,5,10,20,50,100,500,1000]:
    integrationVal = 0.0
    points,weights = np.polynomial.legendre.leggauss(N)
    for i in range(N):
        x = points[i]
        w = weights[i]
        b = getB(x,xMin,xMax,bMin,bMax)
        integrationThing = (bMax-bMin)/(xMax-xMin)
        rhoVal = getVal(xs,x,rho)
        integrationVal += (w*rhoVal/(b*np.tanh(b*0.5))*np.cos(b*t)*integrationThing)
    integrationVals.append(integrationVal)
plt.plot(integrationVals)
plt.plot([x1 for i in range(len(integrationVals))])
plt.show()





















