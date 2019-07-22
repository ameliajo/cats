from phononDists import *
import numpy as np
import matplotlib.pyplot as plt
from math import pi
#import matplotlib.colors as colors
#import matplotlib.cm as cmx

def getP(beta,betas,rho):
    half = int(len(rho)*0.5+1)
    return rho[half]/(betas[half]*betas[half]) if abs(beta) < betas[half] else \
           np.interp(beta,betas,rho)/(2.0*beta*np.sinh(beta*0.5)) 


original_rho   = rho[::-1] + rho[1:]
original_betas = [-x for x in betas][::-1] + betas[1:]
original_P = [getP(beta,original_betas,original_rho) for beta in original_betas]

thin_betas = list(np.linspace(original_betas[0],original_betas[-1],1000))
thin_P = [getP(beta,original_betas,original_rho) for beta in thin_betas]

#plt.plot(original_betas,original_P,linewidth=2)
#plt.plot(thin_betas,thin_P,linewidth=1)
#plt.title('F integrand = P(b)*exp(-b/2)*cos(b*t), t = 62.85 = 2pi/deltaBeta')
#plt.xlabel('beta')

original_integrals = []
thin_integrals = []

t = 62.84985568390593
t_vec = list(np.linspace(0.0,20*t+2,1000))
for i,t in enumerate(t_vec):
    if (i%500==0):
        print(100*i/len(t_vec))
    original_f = [original_P[i]*np.exp(-original_betas[i]*0.5)*np.cos(-original_betas[i]*t) for i in range(len(original_betas))]
    thin_f = [thin_P[i]*np.exp(-thin_betas[i]*0.5)*np.cos(-thin_betas[i]*t) for i in range(len(thin_betas))]

    #plt.plot(original_betas,original_f,label='original beta grid')
    #plt.plot(thin_betas,thin_f,label='thin beta grid')
    original_integrals.append(np.trapz(original_f,x=original_betas))
    thin_integrals.append(np.trapz(thin_f,x=thin_betas))






plt.title('F(t)')
plt.xlabel('t')
plt.plot(t_vec,original_integrals,label='original grid')
plt.plot(t_vec,thin_integrals,label='thin grid')
plt.legend(loc='best')
plt.show()







