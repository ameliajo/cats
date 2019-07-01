import numpy as np
import matplotlib.pyplot as plt
from testRho import *

xMin=-1.0
xMax= 1.0

def getX(b,xMin,xMax,bMin,bMax):
    return (b-bMin)/(bMax-bMin) * (xMax-xMin) + xMin

def getB(x,xMin,xMax,bMin,bMax):
    return (x-xMin)/(xMax-xMin) * (bMax-bMin) + bMin


kbT = 0.025

bs = [x/kbT for x in energies[:]]
bMin, bMax = 0.0, bs[-1]
xs = [getX(b,xMin,xMax,bMin,bMax) for b in bs]

t = 1.0


integrand = []
for i,b in enumerate(bs):
    if b > 0.0:
        integrand.append(rho[i]*np.sin(b*t)/b)
    else:
        integrand.append(rho[i])

print("\nintegrating b1->bmax db'")
print(np.trapz(integrand,x=bs))



integrand = []
for i,x in enumerate(xs):
    b = getB(x,xMin,xMax,bMin,bMax)
    integrationThing = (bMax-bMin)/(xMax-xMin)
    if b > 0.0:
        integrand.append(rho[i]*np.sin(b*t)/b*integrationThing)
    else:
        integrand.append(rho[i]/integrationThing)
print("\nintegrating -1->1 dx")
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
        if b > 1e-6:
            integrationVal += (w*rhoVal/b*np.sin(b*t)*integrationThing)
        else:
            integrationVal += (w*rhoVal*integrationThing)
    integrationVals.append(integrationVal)
plt.plot(integrationVals)
plt.plot([x1 for i in range(len(integrationVals))])
plt.show()
print("\nintegrating with gauss-legendre")
print(integrationVals[-1])
print()



















