import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path.append('../../phononDistributions')
from scipy.interpolate import interp1d
from math import exp
from getSAB_phononExpansion import * 
from getSAB_T1_approx import * 
from colors import *
#from waterDataContinuous import X as rho_x, Q as rho_y; name = 'H in H2O'
from beoData             import X as rho_x, Q as rho_y; name = 'Be in BeO'

rho_x, rho_y = [0.0]+list(rho_x), [0.0]+list(rho_y)
f = interp1d(rho_x,rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
uniform_x = np.linspace(0,rho_x[-1],1*len(rho_x))

betas = np.linspace(0,5,101)
delta = (uniform_x[1]-uniform_x[0])/0.0255

nphon = 10
for a,alpha in enumerate([0.001,0.01,0.1,1.0,5.0]):
    print(alpha)
    sab_nonsym_neg_side = getSAB_phononExpansion(nphon,delta,f(uniform_x),[alpha],betas)
    sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
    plt.plot(betas,sab_nonsym_pos_side, color=misccolors[a], \
             label='alpha = '+str(alpha), linestyle='solid')

    sab = getSAB_T1_approx(delta,f(uniform_x),[alpha],betas)
    plt.plot(betas,[sab[b]*exp(-betas[b]) for b in range(len(betas))],\
             color=misccolors[a],linestyle='dashed', markersize=4)


plt.title('S(a,b) for '+name+'(solid = 10 phonon terms, dashed = approx)')
plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
plt.yscale('log')
plt.show()

