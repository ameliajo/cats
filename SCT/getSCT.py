import matplotlib.pyplot as plt
import numpy as np 
from numpy import sin,cos,sinh,cosh,exp
from math import pi



def getEffectiveTemp(rhoBetas,rho,T):
    P = [rho[b]/(2.0*rhoBetas[b]*sinh(rhoBetas[b]*0.5)) for b in range(len(rho))]
    integrand1 = [P[b]*rhoBetas[b]**2*exp(-rhoBetas[b]) for b in range(len(rho))]
    integrand2 = [P[b]*rhoBetas[b]**2*exp( rhoBetas[b]) for b in range(len(rho))]
    return (np.trapz(integrand1,rhoBetas)+np.trapz(integrand2,rhoBetas))*T*0.5

def SCT(alpha,beta,T,T_eff):
    return (4.0*pi*alpha*T_eff/T)**0.5             * \
            exp(-(alpha-beta)**2/(alpha*T_eff/T)) * \
            exp(-beta)



if __name__=="__main__":
    import sys
    sys.path.append("../phononExpansion");     
    sys.path.append("../phononExpansion/help");
    from getSAB_phononExpansion import *


    sys.path.append("../phononDistributions");     
    #from waterData import rho_f, X as rho_x, Q as rho_y, title
    from beoData   import rho_f, X as rho_x, Q as rho_y, title

    T = 0.0255
    invT = 1.0/T
    rhoBetas = [rhoX_val*invT for rhoX_val in rho_x]

    alphas = [0.01, 10, 100, 300, 500]
    betas  = list(np.linspace(0,100,201))

    T_eff = getEffectiveTemp(rhoBetas,rho_y,T)
    for a,alpha in enumerate(alphas):
        sct_sab = [SCT(alpha,beta,T,T_eff) for beta in betas]
        plt.plot(betas,sct_sab,linestyle='solid')

    plt.title('S(a,b) usign SCT Approximation')
    plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
    plt.yscale('log')
    plt.show()




