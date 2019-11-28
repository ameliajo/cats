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





def simpleGASKET(rhoBetas,rho,time,alphas,betas,flagUnstable = False):
    print('getting F and H...')
    invArea = 1.0/np.trapz(rho,rhoBetas); rho = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)

    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0])/pi for alpha in alphas] # H[0] = debye waller 

    for a,alpha in enumerate(alphas):
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]

        oscBegin = None
        treatOsc   = None

        for b,beta in enumerate(betas):
            i = b+a*len(betas)
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[i] = np.trapz(integrand,time) * alpha_exp[a] 

            if i > 5 and oscBegin == None: 
                if largeDifferences(sab[i-3:i+1]) and slopesDiffer(sab[i-3:i+1]):
                    oscBegin = betas[b]
                    treatOsc   = betas[b-5] 

    return sab,oscBegin,treatOsc





if __name__=="__main__":
    import sys
    sys.path.append("../../phononDistributions"); from waterData import *
    sys.path.append("../../phononExpansion"); from getSAB_phononExpansion import *
    from colors import misccolors

    alphas = [0.001]
    alphas = [0.001]
    betas  = list(np.linspace(0,13,81))
    nbeta = len(betas)
    time = np.linspace(0,96.2,8e2)
    delta = uniform_x[1]-uniform_x[0]

    sab,oscBegin,treatOsc = simpleGASKET(uniform_x[1:],uniform_y[1:],time,alphas,betas)
    print("Just finished gasket")

    for a in range(len(alphas)):
        plt.plot(betas,[sab[b+a*nbeta] for b in range(nbeta)],'ro',label='GASKET',linestyle='solid')
        if oscBegin != None:
            plt.plot([oscBegin,oscBegin],[1e-8,1e-2],'y',label='First detectable osc.')
            plt.plot([treatOsc,treatOsc],[1e-8,1e-2],'g',label='Start treting osc.')

    nphon = 10
    for a in range(len(alphas)):
        sab_neg_side = getSAB_phononExpansion(nphon,delta,uniform_y,[alphas[a]],betas)
        sab_pos_side = [sab_neg_side[b]*exp(-betas[b]) for b in range(nbeta)]
        print("Just finished phonon distribution")
        plt.plot(betas,sab_pos_side,'bo',label='Phonon Expansion')

    plt.yscale('log')
    plt.legend(loc='best')
    plt.show()


