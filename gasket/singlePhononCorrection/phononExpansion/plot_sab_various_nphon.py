import matplotlib.pyplot as plt
import sys
import numpy as np
from math import exp
from contin import contin
sys.path.append('../../phononDistributions')
from waterDataContinuous import X as rho_x
from waterDataContinuous import Q as rho_y
from scipy.interpolate import interp1d

misccolors = ['#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a',\
    '#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2',\
    '#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5']

rho_x = [0.0]+rho_x
rho_y = [0.0]+rho_y
f = interp1d(rho_x,rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
uniform_x = np.linspace(0,rho_x[-1],2*len(rho_x))
uniform_y = f(uniform_x)

alpha = 0.01
alpha = 0.10
alpha = 1.00
betas = np.linspace(0,10,500)
uniform_x = [x/0.0255 for x in uniform_x]


for n,nphon in enumerate([1,2,5,20]):#,4,5,20,30]:#,50]:
    print(nphon)
    sab = contin(nphon,uniform_x[1]-uniform_x[0],uniform_y,[alpha],betas)
    plt.plot(betas,[sab[b]*exp(-betas[b]) for b in range(len(betas))],label='nphon ='+str(nphon),color=misccolors[n])
    plt.plot([-b for b in betas],[sab[b] for b in range(len(betas))],color=misccolors[n])
plt.title('S(a,b) for H in H2O for alpha = '+str(alpha))
plt.xlabel('beta')
plt.ylabel('S(a,b)')
plt.yscale('log')
plt.legend(loc='best')
plt.show()

