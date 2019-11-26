import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../phononDistributions')
from waterDataContinuous import X as rho_x
from waterDataContinuous import Q as rho_y
from scipy.interpolate import interp1d
from start  import *
from convol import *

misccolors = ['#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a',\
    '#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2',\
    '#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5']

def plot_Tn(delta,tn_neg,color,label):
    plt.plot([-delta*i for i in range(len(tn_neg))],tn_neg,color=color,label=label)
    tn_pos = [tn_neg[i]*exp(-delta*i) for i in range(len(tn_neg))]
    plt.plot([ delta*i for i in range(len(tn_pos))],tn_pos,color=color)

f = interp1d([0.0]+rho_x,[0.0]+rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
uniform_x = np.linspace(0,rho_x[-1],3*len(rho_x))
uniform_y = f(uniform_x)
uniform_x = [x/0.0255 for x in uniform_x]
delta = uniform_x[1]-uniform_x[0]

lambda_s, t1 = start(uniform_x,uniform_y,1.0)
tLast = t1[:]; nLast = len(t1)
plot_Tn(delta,t1,misccolors[0],'t1')
nphon = 10 

for i in range(nphon-1):
    nNext = len(t1)+nLast-1
    tNext = convol(t1,tLast,delta,nNext)
    plot_Tn(delta,tNext,misccolors[i+1],'t'+str(i+2))
    tLast = tNext[:]; nLast = nNext


plt.legend(loc='best')
plt.xlim([-2.2*delta*len(t1),2.2*delta*len(t1)])
#plt.ylim([1e-20,1e1])
plt.ylim([-0.1,3.5])
plt.title("Shape of different Tn values for Phonon Expansion Sum")
plt.xlabel("beta'"); plt.ylabel("Tn Values")
#plt.yscale('log');   plt.show()
plt.show()


