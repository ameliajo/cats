import matplotlib.pyplot as plt
import sys
import numpy as np
from math import exp
from contin import contin
sys.path.append('../../phononDistributions')
from waterDataContinuous import X as rho_x
from waterDataContinuous import Q as rho_y
from scipy.interpolate import interp1d


f = interp1d(rho_x,rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
uniform_x = np.linspace(0,rho_x[-1],5*len(rho_x))
uniform_y = f(uniform_x)



alpha = 0.01
betas = np.linspace(0,10,500)
uniform_x = [x/0.0255 for x in uniform_x]
#plt.plot(rho_x,rho_y)
#plt.plot(uniform_x,uniform_y)
#plt.show()


for nphon in [1,2,3,4,5,20]:
    sab = contin(nphon,uniform_x[1]-uniform_x[0],1.0,uniform_y,[alpha],betas)
    plt.plot(betas,[sab[b]*exp(betas[b]) for b in range(len(betas))],label='nphon ='+str(nphon))
plt.title('S(a,b) for H in H2O for alpha = '+str(alpha))
plt.xlabel('beta')
plt.ylabel('S(a,b)')
plt.yscale('log')
plt.legend(loc='best')
plt.show()

