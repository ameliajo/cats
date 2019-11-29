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


def largeDifference(sab1, sab2):
    return abs(sab2-sab1) > (abs(0.5*(sab2+sab1)))

def largeDifferences(sabChunk):
    for i in range(len(sabChunk)-1):
        if not largeDifference(sabChunk[i],sabChunk[i+1]):
            return False
    return True

def slopesDiffer(sabChunk):
    signs = [None]*(len(sabChunk)-1)
    for i in range(1,len(sabChunk)-1):
        if np.sign(sabChunk[i  ]-sabChunk[i-1]) == \
           np.sign(sabChunk[i+1]-sabChunk[i  ]):
            return False
    return True





def simpleGASKET(rhoBetas,rho,time,alphas,betas,useCorrection = False):
    print('getting F and H...')
    invArea = 1.0/np.trapz(rho,rhoBetas); rho = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)

    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0])/pi for alpha in alphas] # H[0] = debye waller 
    unstable = False

    for a,alpha in enumerate(alphas):
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]

        oscBegin = None
        treatOsc = None

        for b,beta in enumerate(betas):
            i = b+a*len(betas)
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[i] = np.trapz(integrand,time) * alpha_exp[a] 

            if useCorrection:
                if i > 10 and oscBegin == None: 
                    if largeDifferences(sab[i-2:i+1]) and slopesDiffer(sab[i-2:i+1]):
                        oscBegin = betas[b]
                        treatOsc = betas[b-10] 
                        unstable = True
                        for j in range(-10,0):
                            sab[b+j] = getSAB_T1_approx(rhoBetas,rho,[alpha],[betas[b+j]])[0]*exp(-betas[b+j])
                if i > 20 and oscBegin == None:
                    if sab[i] < 0:
                        oscBegin = betas[b]
                        treatOsc = betas[b-20]
                        unstable = True
                        for j in range(-20,0):
                            sab[b+j] = getSAB_T1_approx(rhoBetas,rho,[alpha],[betas[b+j]])[0]*exp(-betas[b+j])


                if unstable == True:
                    sab[i] = getSAB_T1_approx(rhoBetas,rho,[alpha],[beta])[0]*exp(-betas[b])

    return sab,H,F,oscBegin,treatOsc





if __name__=="__main__":
    import sys
    sys.path.append("../../phononDistributions"); from waterData import *
    sys.path.append("../../phononExpansion"); from getSAB_phononExpansion import *
    from colors import misccolors

    alphas = [0.001]
    betas = list(np.linspace(0,13,201)); nbeta = len(betas)
    time = np.linspace(0,96.2,8e2)
    delta = uniform_x[1]-uniform_x[0]

    #sab_approx = [getSAB_T1_approx(uniform_x[1:],uniform_y[1:],[alphas[0]],[betas[b]])[0]*exp(-betas[b]) for b in range(len(betas))]
    #plt.plot(betas,sab_approx,'yo',markersize=5,linestyle='solid')

    #for i in range(len(betas)):
    #    print(betas[i])
    #print(len(betas))
    #exit()

    sab_no_correction,  H, F, oscBegin, treatOsc = \
              simpleGASKET(uniform_x[1:],uniform_y[1:],time,alphas,betas,False)
    sab_yes_correction, H, F, oscBegin, treatOsc = \
              simpleGASKET(uniform_x[1:],uniform_y[1:],time,alphas,betas,True)
    print("Just finished gasket")

    for a in range(len(alphas)):
        plt.plot(betas,[sab_no_correction[b+a*nbeta] for b in range(nbeta)],\
                 'r',label='no correction',linestyle='solid',linewidth=4)
        plt.plot(betas,[sab_yes_correction[b+a*nbeta] for b in range(nbeta)],\
                 'y',label='with correction',linestyle='solid',linewidth=3)

        if oscBegin != None:
            plt.plot([oscBegin,oscBegin],[1e-9,1e-2],'b',label='First detectable osc.')
            plt.plot([treatOsc,treatOsc],[1e-9,1e-2],'g',label='Start treting osc.')

    #nphon = 10
    #for a in range(len(alphas)):
    #    sab_neg_side = getSAB_phononExpansion(nphon,delta,uniform_y,[alphas[a]],betas)
    #    sab_pos_side = [sab_neg_side[b]*exp(-betas[b]) for b in range(nbeta)]
    #    print("Just finished phonon distribution")
    #    plt.plot(betas,sab_pos_side,'k',label='Phonon Expansion')

    plt.yscale('log')
    plt.legend(loc='best')
    plt.show()


