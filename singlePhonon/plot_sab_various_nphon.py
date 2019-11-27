import matplotlib.pyplot as plt
import sys
import numpy as np
from math import exp
sys.path.append('../phononDistributions')
sys.path.append('../phononExpansion')
sys.path.append('../phononExpansion/help')
from getSAB_phononExpansion import *
from colors import misccolors

#from waterData    import rho_f, X as rho_x, title
#from beoData      import rho_f, X as rho_x, title
from graphiteData import rho_f, X as rho_x, title

uniform_x = np.linspace(0,rho_x[-1],2*len(rho_x))
uniform_y = rho_f(uniform_x)

alpha = 0.01
alpha = 0.10
#alpha = 1.00

betas = np.linspace(0,10,500)
uniform_x = [x/0.0255 for x in uniform_x]

for n,nphon in enumerate([1,2,5]):#,4,5,20,30]:#,50]:
    print(nphon)
    sab_nonsym_neg_side = getSAB_phononExpansion(nphon,uniform_x[1]-uniform_x[0],uniform_y,[alpha],betas)
    sab_nonsym_pos_side = [sab_nonsym_neg_side[b]*exp(-betas[b]) for b in range(len(betas))]
    plt.plot(betas,sab_nonsym_pos_side,label='nphon ='+str(nphon),color=misccolors[n])
    plt.plot([-b for b in betas],sab_nonsym_neg_side,color=misccolors[n])

plt.title('Effect of nphon on S(a,b) for '+title+' for alpha = '+str(alpha))
plt.xlabel('beta');     plt.ylabel('S(a,b)'); plt.yscale('log')
plt.legend(loc='best'); plt.show()

