import numpy as np
from getF import *
from getG import *
from getLambda_s import *
import matplotlib.pyplot as plt
from math import pi

kbT = 296.0*8.6173332e-5
betas = [x/kbT for x in energies[:]]
invArea = 1.0/np.trapz(rho,x=betas)
rho = [invArea*x for x in rho]

lambda_s = integrate_Lambda_XsGaussLegendre(betas[1:],rho,100)+getSmallContrib(betas,rho)
print('\n',lambda_s,'\n')

alpha = 0.1
beta = 0.2
alpha = 0.01
beta = 0.006



def plot_integrand_and_integration(alpha,beta,betas,rho):
    tVec = np.linspace(0,100,300)

    contribs, integrations = [], []
    for i,t in enumerate(tVec):
        F, G = getF(betas,rho,t,100), getF(betas,rho,t,100)
        contribs.append( (np.exp(-alpha*lambda_s)/pi) * np.exp(alpha*F) * np.cos(beta*t - alpha*G) )
        #if i > 5: integrations.append(np.trapz(contribs,x=tVec[:i+1]))
        #else:     integrations.append(0.0)
    plt.plot(tVec,contribs)
    #plt.plot(tVec,integrations)
    plt.show()
 
#plot_integrand_and_integration(alpha,beta,betas,rho)

def integrate_GaussLaguerre(betas,rho,N):
    contribs = 0.0
    points,weights = np.polynomial.laguerre.laggauss(N)
    for i in range(N):
        t = points[i]
        w = weights[i]

        F, G = getF(betas,rho,t,50), getF(betas,rho,t,10)

        contribs += w * (np.exp(-alpha*lambda_s)/pi) * np.exp(alpha*F) * np.cos(beta*t - alpha*G)
    return contribs


#for N in [2,5,10,20,50,100]:
#    print(integrate_GaussLaguerre(betas,rho,N))

#plot_integrand_and_integration(alpha,beta,betas,rho)












