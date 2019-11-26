import sys; sys.path.append('./help/')
import numpy as np
from getT1       import *
from convolution import *
from interpolate import *
from math        import exp, factorial

def getSAB_phononExpansion(nphon, delta, rho, alphas, betas, wgt = 1.0):
    betaGrid     = [delta*i for i in range(len(rho))]
    lambda_s, t1 = getT1(betaGrid,rho,wgt) # This is T1(-beta) 
    sab = [0.0]*len(alphas)*len(betas)
    xa  = [1.0]*len(alphas)
    tnow  = [0.0]*(len(t1)*nphon)
    tlast = [0.0]*(len(t1)*nphon)
    nLast = len(t1)
    for i in range(len(t1)):
        tnow[i]  = t1[i]
        tlast[i] = t1[i]
    lambda_alpha = [lambda_s*alpha for alpha in alphas]
    exp_lambda_alpha = [exp(-lambda_alpha_val) for lambda_alpha_val in lambda_alpha]

    for n in range(nphon):
        if n > 0:
            # Convol takes T1(-beta) and T2(-beta) and gives you T3(-beta)
            nNext = len(t1)+nLast-1
            tnow = convol(t1, tlast, delta, nNext)
        inv_n = 1.0/(n+1)
        for a in range(len(alphas)):
            xa[a] *= lambda_alpha[a] * inv_n
            exx    = exp_lambda_alpha[a]*xa[a]

            for b in range(len(betas)):
                sab[b+a*len(betas)] += exx * interpolate(tnow,betas[b],delta)
                #add = exp(-alphas[a]*lambda_s) * \
                #      interpolate(tnow,betas[b],delta) * \
                #      (alphas[a]*lambda_s)**(n+1) / factorial(n+1)
                #if add > 1e-30: sab[b+a*len(betas)] += add
        if n == 0:
            continue
        for i in range(nNext):
            tlast[i] = tnow[i]
        nLast = nNext
    return sab

if __name__=='__main__':
    sys.path.append('../../phononDistributions')
    from waterDataContinuous import X as rho_x
    from waterDataContinuous import Q as rho_y
    from colors              import misccolors
    from scipy.interpolate   import interp1d
    import matplotlib.pyplot as plt

    rho_x, rho_y = [0.0]+list(rho_x), [0.0]+list(rho_y)
    f = interp1d(rho_x,rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
    uniform_x = np.linspace(0,rho_x[-1],1*len(rho_x))
    delta = (uniform_x[1]-uniform_x[0])/0.0255
    betas = np.linspace(0,5,101)

    nphon = 10
    for a,alpha in enumerate([0.001,0.01,0.1,1.0]):
        print('alpha = '+str(alpha))
        sab_nonsym_neg_side = getSAB_phononExpansion(nphon,delta,f(uniform_x),[alpha],betas)
        sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
        plt.plot(betas,sab_nonsym_pos_side,color=misccolors[a],\
                 label='alpha = '+str(alpha),linestyle='solid')
    plt.title('S(a,b) for H in H2O, for full phonon expansion with '+str(nphon)+' terms')
    plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
    plt.yscale('log')
    plt.show()






