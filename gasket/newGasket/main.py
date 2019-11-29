import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np 
from numpy import sin,cos,sinh,cosh,exp
from math import pi
from get_F_H import *



def simpleGASKET(rhoBetas,rho,time,alphas,betas,flagUnstable = False):
    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho     = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)

    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0])/pi for alpha in alphas] # H[0] = debye waller 
    beta_exp  = [exp(beta*0.5)       for beta  in betas ] # This is to turn 
                                         # the scattering law to be symmetric

    #unstable = False
    for a,alpha in enumerate(alphas):
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]
        #numNeg = 0
        for b,beta in enumerate(betas):
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[b+a*len(betas)] = np.trapz(integrand,time) * alpha_exp[a] * beta_exp[b]
            #if b+a*len(betas) < 2: continue
            
            #if sab[b+a*len(betas)-2] < 0 or \
            #   sab[b+a*len(betas)-1] < 0 or \
            #   sab[b+a*len(betas)  ] < 0:   \
            #       unstable = True
            #if sab[b+a*len(betas)-2] > 0 and \
            #   sab[b+a*len(betas)-1] > 0 and \
            #   sab[b+a*len(betas)  ] > 0:    \
            #       unstable = False

            #if unstable:
            #    plt.plot(beta,sab[b+a*len(betas)],'ro')
            #else:
            #    plt.plot(beta,sab[b+a*len(betas)],'bo')


    return sab,H,F





if __name__=="__main__":
    import sys
    sys.path.append("../../phononDistributions"); from waterData import *
    from colors import misccolors

    invT = 1.0/0.025
    rhoBetas = [rhoX_val*invT for rhoX_val in X]

    alphas = [0.001,0.1,1.0,5.0,8.0]#,10.0,20.0]
    alphas = [0.001,0.1,1.0]#5.0,8.0]#,10.0,20.0]
    alphas = [0.001]
    #betas  = list(np.linspace(0,10,41))
    betas1  = list(np.linspace(3,9,21))

    NT1 = 5e3; time1 = np.linspace(0,94,NT1)

    sab1,H,F = simpleGASKET(rhoBetas,Q,time1,alphas,betas1)
    plt.yscale('log')
    betas2 = betas1 + [4.1,4.13,4.15,4.18]
    betas2.sort()
    sab2,H,F = simpleGASKET(rhoBetas,Q,time1,alphas,betas2)

    for a in range(len(alphas)):
        plt.plot(betas1,[sab1[b+a*len(betas1)] for b in range(len(betas1))],\
                 color=misccolors[a],linewidth=1.5,alpha=0.8)
        plt.plot(betas2,[sab2[b+a*len(betas2)] for b in range(len(betas2))],\
                 color=misccolors[a],linewidth=1.5,alpha=0.8)



    plt.yscale('log')
    plt.show()


