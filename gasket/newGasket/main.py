from get_F_H import *
import numpy as np
from numpy import sin,cos,sinh,cosh,exp
import matplotlib.pyplot as plt



def gasket(rhoBetas,rho,time,alphas,betas,experiment=False):
    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho = [invArea * x for x in rho]

    F, H = get_F_H(rhoBetas,rho,time,experiment)
    print("got F and H")
    #plt.plot(time,F,label="F")
    #plt.plot(time,H,label="H")
    #plt.legend(loc='best')
    #plt.show()
    lambda_s  = H[0]
    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*lambda_s) for alpha in alphas]
    for a,alpha in enumerate(alphas):
        print(a)
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]
        for b,beta in enumerate(betas):
            #integrand = [exp(alpha*H[t]) * cos(beta*time[t]+alpha*F[t]) for t in range(len(time))]

            #integrand = [exp(alpha*H[t]) * cos(alpha*F[t]) * cos(beta*time[t]) \
            #           - exp(alpha*H[t]) * sin(alpha*F[t]) * sin(beta*time[t]) \
            #           - cos(beta*time[t]) for t in range(len(time))]
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[b+a*len(betas)] = np.trapz(integrand,time)*alpha_exp[a]
    return sab



from waterDataContinuous import *
invT = 1.0/0.025
rhoBetas = [rhoX_val*invT for rhoX_val in X]
rho = Q[:]

colors = ['#35353C', '#4D5E43', '#FFC132', '#FF7A32', '#FF3232']

alphas = [0.001,0.1,1.0,2.0,5.0]
betas  = list(np.linspace(0,10,101))

NT = 8e3
time1   = np.linspace(0,100,NT)
time2 = list(np.linspace(0.0,8.0,int(NT)*0.50)) + \
        list(np.linspace(8.0001,20.0,int(NT)*0.40)) + \
        list(np.linspace(20.0001,100.0,int(NT)*0.10)) 
sab1 = gasket(rhoBetas,rho,time1,alphas,betas)
sab2 = gasket(rhoBetas,rho,time1,alphas,betas,experiment=True)

for a in range(len(alphas)):
    plt.plot(betas,[sab1[b+a*len(betas)] for b in range(len(betas))],\
             linestyle='solid',color=colors[a])
    plt.plot(betas,[sab2[b+a*len(betas)] for b in range(len(betas))],\
            linestyle='dashed',color=colors[a])



plt.yscale('log')
plt.show()







