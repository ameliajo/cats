import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,sinh,cosh,exp
from get_F_H import *
from math import pi


def simpleGASKET(rhoBetas,rho,time,alphas,betas):
    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho     = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)
    print("got F and H")
    #plt.plot(time,F,label="F"); plt.plot(time,H,label="H")
    #plt.legend(loc='best');     plt.show(); exit()
    inv_pi    = 1.0/pi
    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0])*inv_pi for alpha in alphas] # H[0] = debye waller 
    beta_exp  = [exp(beta*0.5)           for beta  in betas ] # This is to turn 
                                         # the scattering law to be symmetric
    for a,alpha in enumerate(alphas):
        print((a+1)/len(alphas)*100)
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]

        numNeg = 0

        for b,beta in enumerate(betas):
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[b+a*len(betas)] = np.trapz(integrand,time) * alpha_exp[a] * beta_exp[b]
            if sab[b+a*len(betas)] < 0: numNeg += 1
            #    print(a,b,sab[b+a*len(betas)])
        print("Found",numNeg,"negative values")

    #print("-----",len(list(filter(lambda x: x < 0, sab))))
    #print(" ")
    return sab,H,F





















if __name__=="__main__":
    from waterDataContinuous import *
    invT = 1.0/0.025
    rhoBetas = [rhoX_val*invT for rhoX_val in X]
    rho = Q[:]

    colors = ['#35353C', '#4D5E43', '#FFC132', '#FF7A32', '#FF3232']
    
    alphas = [0.001,0.1,1.0,2.0,5.0]
    betas  = list(np.linspace(0,8,81))

    NT = 5e3
    time1   = np.linspace(0,100,NT)
    #time2 = list(np.linspace(0.0,8.0,int(NT)*0.50)) + \
    #        list(np.linspace(8.0001,20.0,int(NT)*0.40)) + \
    #        list(np.linspace(20.0001,100.0,int(NT)*0.10)) 
    sab1,H,F = simpleGASKET(rhoBetas,rho,time1,alphas,betas)
    #sab2,H,F = simpleGASKET(rhoBetas,rho,time2,alphas,betas)

    for a in range(len(alphas)):
        plt.plot(betas,[sab1[b+a*len(betas)] for b in range(len(betas))],\
                 linestyle='solid',color=colors[a])
        #plt.plot(betas,[sab2[b+a*len(betas)] for b in range(len(betas))],\
        #        linestyle='dashed',color=colors[a])

    plt.yscale('log')
    plt.show()







