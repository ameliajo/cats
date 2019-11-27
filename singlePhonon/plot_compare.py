import numpy as np
import matplotlib.pyplot as plt
import sys 
sys.path.append('../phononDistributions')
sys.path.append('../phononExpansion')
sys.path.append('../phononExpansion/help')
from math import exp
from getSAB_phononExpansion import * 
from getSAB_T1_approx import * 
from colors import misccolors

#from beoData import rho_f, X as rho_x, title
from waterData import rho_f, X as rho_x, title
#from graphiteData import rho_f, X as rho_x, title

uniform_x = np.linspace(0,rho_x[-1],1*len(rho_x))
uniform_y = rho_f(uniform_x)
delta = (uniform_x[1]-uniform_x[0])/0.0255
uniform_x = [x/0.0255 for x in uniform_x]

nphon = 10
betas = np.linspace(0,25,501)
for a,alpha in enumerate([0.001]):#,0.01,0.1,1.0]):
    print(alpha)
    sab_nonsym_neg_side = getSAB_phononExpansion(nphon,delta,uniform_y,[alpha],betas)
    sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
    plt.plot(betas,sab_nonsym_pos_side, color=misccolors[a], \
             label='alpha = '+str(alpha), linestyle='solid')
    #print(sab_nonsym_neg_side[:5])
    #print(sab_nonsym_pos_side[:5])

    sab = getSAB_T1_approx(uniform_x,uniform_y,[alpha],betas)
    plt.plot(betas,[sab[b]*exp(-betas[b]) for b in range(len(betas))],\
             color=misccolors[a],linestyle='dashed', markersize=4)


plt.title('S(a,b) for '+title+'(solid = 10 phonon terms, dashed = approx)')
plt.xlabel('beta'); plt.ylabel('S(a,b)'); plt.legend(loc='best')
plt.yscale('log')
plt.show()

