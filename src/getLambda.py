import numpy as np
from testRho import *
import matplotlib.pyplot as plt

kbT = 0.025

fullBeta = [E/kbT for E in fullE]
invArea = 1.0/np.trapz(fullRho, x=fullBeta)
fullRho = [x*invArea*2 for x in fullRho]
rho = fullRho[len(rho)-1:]


P = [fullRho[b]/(2*fullBeta[b]*np.sinh(fullBeta[b]*0.5)) if fullBeta[b] != 0.0 else fullRho[b]/(spacing/kbT)**2 for b in range(len(fullE))]


integrand = [P[b]*np.exp(-fullBeta[b]*0.5) for b in range(len(fullE))]
lambda_s = np.trapz(integrand,x=fullBeta)

betas = [E/kbT for E in energies]
integrand = [2.0*rho[1]/(spacing/kbT)**2]
for b in range(1,len(energies)):
    integrand.append(rho[b]/(betas[b]) * (np.cosh(betas[b]*0.5)/np.sinh(betas[b]*0.5))*np.exp(betas[b]))

lambda_s = np.trapz([integrand[b]*np.exp(-betas[b]) for b in range(len(betas))],x=betas)



plt.plot([lambda_s for x in range(100)])
lambda_list = []
for N in range(1,100):
    points,weights = np.polynomial.laguerre.laggauss(N)
    lambda_s = 0.0
    for i in range(N):
        lambda_s += weights[i]*getVal(betas,points[i],integrand)
    lambda_list.append(lambda_s)
plt.plot(lambda_list)
plt.show()












