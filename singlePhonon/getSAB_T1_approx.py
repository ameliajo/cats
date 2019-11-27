import sys; sys.path.append('../phononExpansion/help/')
import numpy as np
from getT1       import *
from convolution import *
from interpolate import *
from math        import exp, factorial


def getSAB_T1_approx(betaGrid, rho, alphas, betas):
    lambda_s, t1 = getT1(betaGrid,rho,1.0)
    sab  = [0.0]*len(alphas)*len(betas)
    lambda_alpha = [lambda_s*alpha for alpha in alphas]
    exp_lambda_alpha = [exp(-lambda_alpha_val) for lambda_alpha_val in lambda_alpha]

    for a in range(len(alphas)):
        for b in range(len(betas)):
            sab[b+a*len(betas)] = interpolate(betaGrid,t1,betas[b]) * \
                                  (1.0 - exp_lambda_alpha[a])
    return sab


if __name__=='__main__':
    sys.path.append('../phononDistributions')
    from waterData import X as rho_x
    from waterData import Q as rho_y
    from colors              import misccolors
    from scipy.interpolate   import interp1d
    import matplotlib.pyplot as plt

    rho_x, rho_y = [0.0]+list(rho_x), [0.0]+list(rho_y)
    f = interp1d(rho_x,rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
    uniform_x = np.linspace(0,rho_x[-1],1*len(rho_x))
    delta = (uniform_x[1]-uniform_x[0])/0.0255
    betas = np.linspace(0,25,101)

    betaGrid     = [delta*i for i in range(len(rho_x))]

    for a,alpha in enumerate([0.001,0.01,0.1,1.0]):
        print('alpha = '+str(alpha))
        sab_nonsym_neg_side = getSAB_T1_approx(betaGrid,f(uniform_x),[alpha],betas)
        sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
        plt.plot(betas,sab_nonsym_pos_side,color=misccolors[a],\
                 label='alpha = '+str(alpha),linestyle='solid')
    plt.title('S(a,b) for H in H2O, for T1 Approximation')
    plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
    plt.yscale('log')
    plt.show()






