from getDebyeWaller import getDebyeWaller
from normalize      import normalize 
from math import sinh, exp



def start(betas, rho, wgt):
    # Returns Debye-Waller factor lambda_s and the first value in the phonon
    # expansion, T1. Note that this is T1(-beta)
    P = [rho[1]/betas[1]**2] + \
        [rho[b]/(2.0*betas[b]*sinh(betas[b]*0.5)) for b in range(1,len(rho))]
    P = normalize(betas,P,wgt)
    lambda_s = getDebyeWaller(betas,P)
    T1 = [P[b]*exp(betas[b]*0.5)/lambda_s for b in range(len(betas))]
    return lambda_s, T1



if __name__=='__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    import sys; sys.path.append('../../phononDistributions')
    from waterDataContinuous import X as rho_x, Q as rho_y
    from scipy.interpolate import interp1d

    misccolors = ['#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a',\
        '#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2',\
        '#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5']

    def plot_Tn(delta,tn_neg,color,label):
        plt.plot([-delta*i for i in range(len(tn_neg))],tn_neg,color=color,label=label)
        tn_pos = [tn_neg[i]*exp(-delta*i) for i in range(len(tn_neg))]
        plt.plot([ delta*i for i in range(len(tn_pos))],tn_pos,color=color)

    f = interp1d([0.0]+rho_x,[0.0]+rho_y,bounds_error=False,fill_value=0.0,kind='cubic')
    uniform_x = np.linspace(0,rho_x[-1]+1e-5,3*len(rho_x))
    uniform_y = f(uniform_x)
    uniform_x = [x/0.0255 for x in uniform_x]
    delta = uniform_x[1]-uniform_x[0]
    lambda_s, t1 = start(uniform_x,uniform_y,1.0)

