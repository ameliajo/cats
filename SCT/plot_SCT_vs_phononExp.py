import matplotlib.pyplot as plt
import sys
from getSCT import *
sys.path.append("../phononDistributions"); 
sys.path.append("../phononExpansion");     
sys.path.append("../phononExpansion/help");
from getSAB_phononExpansion import *
from colors import misccolors


#from waterData    import rho_f, X as rho_x, Q as rho_y, title
from beoData    import rho_f, X as rho_x, Q as rho_y, title


if __name__=="__main__":

    T = 0.0255
    invT = 1.0/T
    rhoBetas = [rhoX_val*invT for rhoX_val in rho_x]

    alphas = [0.01]#, 10, 100]
    betas  = list(np.linspace(0,50,51))

    T_eff = getEffectiveTemp(rhoBetas,rho_y,T)
    for a,alpha in enumerate(alphas):
        sct_sab = [SCT(alpha,beta,T,T_eff) for beta in betas]
        plt.plot(betas,sct_sab,linestyle='dashed',color=misccolors[a])

    uniform_x = np.linspace(0,rho_x[-1],1*len(rho_x))
    delta = (uniform_x[1]-uniform_x[0])/0.0255
    nphon = 30
    for a,alpha in enumerate(alphas):
        print('alpha = '+str(alpha))
        sab_nonsym_neg_side = getSAB_phononExpansion(nphon,delta,rho_f(uniform_x),[alpha],betas)
        sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
        plt.plot(betas,sab_nonsym_pos_side,color=misccolors[a],\
                 label='alpha = '+str(alpha),linestyle='solid')

    plt.title('S(a,b) for '+title+' comparing phonon sum (solid) with SCT (dashed)')
    plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
    plt.yscale('log')
    plt.show()




