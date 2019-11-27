import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np 
from numpy import sin,cos,sinh,cosh,exp
from math import pi
import sys
sys.path.append('../newGasket/')
sys.path.append('../../singlePhonon/')
sys.path.append('../../phononExpansion/')
sys.path.append('../../phononExpansion/help')
from getSAB_T1_approx import *
from get_F_H import *

def areWeUnstable(a,b,betas,sab,alpha_H_exp,cos_alpha_F,sin_alpha_F,time,alpha_exp):
    F1_x = betas[b-2]; F1_y = sab[b+a*len(betas)-2]
    F2_x = betas[b-1]; F2_y = sab[b+a*len(betas)-1]
    F3_x = betas[b  ]; F3_y = sab[b+a*len(betas)  ]

    if F2_y > F1_y and F3_y < F2_y: # If we're not monotonically increasing or decreasing
        N = 20
        sab_i_vec  = []
        num_osc = 0
        for i in range(N):
            beta_i = i*(F3_x-F2_x)/N + F2_x
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta_i*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta_i*time[t]) \
                       - cos(beta_i*time[t]) for t in range(len(time))]

            sab_at_beta_i = np.trapz(integrand,time) * alpha_exp[a] #* exp(beta_i*0.5)
            sab_i_vec.append(sab_at_beta_i)
            if len(sab_i_vec) > 2:
                if sab_i_vec[-2] - sab_i_vec[-1] > 0.6 * (F3_y-F2_y) and \
                   sab_i_vec[-2] > sab_i_vec[-3]:
                    num_osc += 1
            if sab_at_beta_i < 0:
                return True
        if num_osc > 2:
            return True
    return False




def simpleGASKET(rhoBetas,rho,time,alphas,betas,flagUnstable = False):
    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho     = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)
    #plt.plot(time,F); plt.plot(time,H)
    #plt.show(); exit()

    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0])/pi for alpha in alphas] # H[0] = debye waller 
    #beta_exp  = [exp(beta*0.5)       for beta  in betas ] # This is to turn 
                                         # the scattering law to be symmetric

    unstable = False
    for a,alpha in enumerate(alphas):
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]
        numNeg = 0
        for b,beta in enumerate(betas):
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[b+a*len(betas)] = np.trapz(integrand,time) * alpha_exp[a] 

            if b+a*len(betas) < 2: continue
            

            if not unstable:
                unstable = areWeUnstable(a,b,betas,sab,alpha_H_exp,cos_alpha_F,sin_alpha_F,time,alpha_exp)
                if unstable == True:
                    for i in range(-5,0):
                        sab_approx = getSAB_T1_approx(rhoBetas,rho,[alphas[a]],[betas[b+i]])[0]*exp(-betas[b])
                        plt.plot(betas[b+i],sab_approx,'yo')
 
            if unstable:
                plt.plot(beta,sab[b+a*len(betas)],'ro')
                sab_approx = getSAB_T1_approx(rhoBetas,rho,[alphas[a]],[betas[b]])[0]*exp(-betas[b])
                plt.plot(beta,sab_approx,'yo')
            else:
                plt.plot(beta,sab[b+a*len(betas)],'bo')


    return sab,H,F





if __name__=="__main__":
    import sys
    sys.path.append("../../phononDistributions"); from waterData import *
    sys.path.append("../../phononExpansion"); from getSAB_phononExpansion import *
    from colors import misccolors



    kbT = 8.61733e-5 * 293.0
    invT = 1.0/kbT
    #rhoBetas = [rhoX_val*invT for rhoX_val in X]

    alphas = [0.001]
    betas  = list(np.linspace(0,5,51))

    NT1 = 1e4; time1 = np.linspace(0,120,NT1)

    sab,H,F = simpleGASKET(uniform_x[1:],uniform_y[1:],time1,alphas,betas)
    print("Just finished gasket")

    for a in range(len(alphas)):
        plt.plot(betas,[sab[b+a*len(betas)] for b in range(len(betas))],\
                 color=misccolors[a],linewidth=1.5,alpha=0.8)


        sab_nonsym_neg_side = getSAB_phononExpansion(5,uniform_x[1]-uniform_x[0],uniform_y,[alphas[a]],betas)
        print("Just finished phonon distribution")
        sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
        plt.plot(betas,sab_nonsym_pos_side,'go')

    plt.yscale('log')
    plt.show()


