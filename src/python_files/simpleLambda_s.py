import numpy as np

def getVal(energies,val,vec):
    for i in range(len(energies)-1):
        if energies[i] <= val and val <= energies[i+1]:
            m = (vec[i+1]-vec[i])/(energies[i+1]-energies[i])
            b = vec[i]
            return m*(val-energies[i])+b
    return 0.0
 
rho = [0, 0.09, 0.16, 0.21, 0.24, 0.25, 0.24, 0.21, 0.16, 0.09, 0]
betas = [0.1*i for i in range(len(rho))]
invArea = 1.0/np.trapz(rho,x=betas)
#print("area",1.0/invArea)
rho = [invArea*x for x in rho]

P = [rho[1]/(betas[1]*betas[1])] + [rho[i]/(2.0*betas[i]*np.sinh(betas[i]*0.5)) for i in range(1,len(rho))]
#print(P)
integrand = [P[i]*2*np.cosh(betas[i]*0.5) for i in range(len(rho))]
print(np.trapz(integrand,x=betas))


fullBetas = [-b for b in betas][::-1] + betas[1:]
print(fullBetas)
xs = [b/betas[-1] for b in fullBetas]
fullRho = rho[::-1] + rho[1:]
fullP = P[::-1] + P[1:]
integrationVal = 0.0
N = 10
points,weights = np.polynomial.legendre.leggauss(N)
for i in range(N):
    x = points[i]
    w = weights[i]
    b = betas[-1]*x
    integrationThing = betas[-1]
    pVal = getVal(xs,x,fullP)
    integrationVal += w*pVal*np.exp(-b*0.5)*integrationThing
print(integrationVal)






