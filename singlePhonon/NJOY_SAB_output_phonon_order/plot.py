import matplotlib.pyplot as plt
from BeO_nphon_1   import temp as temp1,   alphas as alphas1,   betas as betas1,   sab as sab1
from BeO_nphon_3   import temp as temp3,   alphas as alphas3,   betas as betas3,   sab as sab3
from BeO_nphon_5   import temp as temp5,   alphas as alphas5,   betas as betas5,   sab as sab5
from BeO_nphon_10  import temp as temp10,  alphas as alphas10,  betas as betas10,  sab as sab10
from BeO_nphon_100 import temp as temp100, alphas as alphas100, betas as betas100, sab as sab100
import sys
sys.path.append('../../../phononDistributions')
from beoData import *

assert(temp1==temp3==temp5==temp10==temp100)
assert(alphas1==alphas3==alphas5==alphas10==alphas100)
assert(betas1==betas3==betas5==betas10==betas100)
temp   = temp1
alphas = alphas1[:]
betas  = betas1[:]



def plot_sab(a,alphas,betas,sab,label,linewidth):
    to_plot_x = [] 
    to_plot_y = [] 
    betaMax = 10
    for b, beta in enumerate(betas):
        if beta < betaMax:
            to_plot_x.append(beta)
            to_plot_y.append(sab[b+a*len(betas)])
        else:
            break
    plt.plot(to_plot_x,to_plot_y,label=label,linewidth=linewidth)


for a in range(len(alphas)):
    betaMax = 10
    plt.figure()
    plt.subplot(211)
    plt.title('BeO alpha = '+str(alphas[a]))
    plt.plot([x/0.0255 for x in X],Q)
    plt.xlim([betas[0],betaMax])

    plt.subplot(212)
    plot_sab(a,alphas,betas,sab1,'n=1',3)
    plot_sab(a,alphas,betas,sab3,'n=3',2.0)
    #plot_sab(a,alphas,betas,sab5,'n=5',1.5)
    #plot_sab(a,alphas,betas,sab10,'n=10',1.0)
    plot_sab(a,alphas,betas,sab100,'n=100',1.0)
    plt.xlim([betas[0],betaMax])
    plt.legend(loc='best')

    plt.yscale('log')
    plt.show()

