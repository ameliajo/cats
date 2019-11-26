import matplotlib.pyplot as plt
import sys; import numpy as np
sys.path.append('../../phononDistributions'); 
sys.path.append('./help/')
from colors      import misccolors
from getT1       import *
from convolution import *
#from beoData import rho_f, X as rho_x, title
#from waterData import rho_f, X as rho_x, title
from graphiteData import rho_f, X as rho_x, title



def plot_Tn(delta,tn_neg,color,label):
    plt.plot([-delta*i for i in range(len(tn_neg))],tn_neg,color=color,label=label)
    tn_pos = [tn_neg[i]*exp(-delta*i) for i in range(len(tn_neg))]
    plt.plot([ delta*i for i in range(len(tn_pos))],tn_pos,color=color)

uniform_x = np.linspace(0,rho_x[-1],3*len(rho_x))
uniform_y = rho_f(uniform_x)
uniform_x = [x/0.0255 for x in uniform_x]
delta = uniform_x[1]-uniform_x[0]

lambda_s, t1 = getT1(uniform_x,uniform_y,1.0)
tLast = t1[:]; nLast = len(t1)
plot_Tn(delta,t1,misccolors[0],'t1')
nphon = 10 

for i in range(nphon-1):
    print('on expansion #'+str(i))
    nNext = len(t1)+nLast-1
    tNext = convol(t1,tLast,delta,nNext)
    plot_Tn(delta,tNext,misccolors[i+1],'t'+str(i+2))
    tLast = tNext[:]; nLast = nNext


plt.legend(loc='best')
plt.title("Various Tn values for Phonon Expansion Sum, ("+title+")")
plt.xlabel("beta'"); plt.ylabel("Tn Values")
plt.yscale('log');   plt.show()


