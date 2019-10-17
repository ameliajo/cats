import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.pyplot import figure
import matplotlib
import numpy as np
from numpy import sin,cos,sinh,cosh,exp
from get_F_H import *
from math import pi

def prepPlot(vector):
    cnorm = colors.Normalize(vmin=0,vmax=len(vector))
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) 
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vector))])
    colorBar = plt.contourf([[0,0],[0,0]], vector, cmap=mymap)
    return scalarMap, colorBar

def prepPlot2(alphas):
    scalarMap, colorBar = prepPlot(alphas); plt.clf()
    colorList = [colors.rgb2hex(scalarMap.to_rgba(i)) for i in range(len(alphas))]
    mymap = colors.ListedColormap(colorList)
    colorBar = plt.contourf([[0,0],[0,0]], list(range(len(alphas)+1)), cmap=mymap)
    plt.clf()
    matplotlib.rcParams.update({'font.size': 15})
    figure(num=None, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    plt.title('H in H2O Scattering Law (Direct Method)')
    plt.ylabel('S(a,b) (Non-symmetric)')
    plt.xlabel('beta')
    plt.yscale('log')

    cbar = plt.colorbar(colorBar)
    cbar.ax.get_yaxis().set_ticks([])
    plt.rcParams.update({'font.size': 14})
    for j, lab in enumerate(["    "+str(a) for a in alphas[:]]):
        cbar.ax.text(1.7, (2 * j + 1) / 13.9, lab, ha='center', va='center')
    cbar.ax.get_yaxis().labelpad = 40
    cbar.ax.set_ylabel('alpha')
    return scalarMap




def simpleGASKET(rhoBetas,rho,time,alphas,betas):
    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho     = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)
    #F[-1], H[-1] = 0.0, 0.0
    print("got F and H")

    inv_pi    = 1.0/pi
    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0])*inv_pi for alpha in alphas] # H[0] = debye waller 
    beta_exp  = [exp(beta*0.5)           for beta  in betas ] # This is to turn 
                                         # the scattering law to be symmetric
    for a,alpha in enumerate(alphas):
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]
        numNeg = 0
        for b,beta in enumerate(betas):
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[b+a*len(betas)] = np.trapz(integrand,time) * alpha_exp[a] * beta_exp[b]
            if sab[b+a*len(betas)] < 0: numNeg += 1
        print("Found",numNeg,"negative values")

    return sab,H,F





















if __name__=="__main__":
    from waterDataContinuous import *
    invT = 1.0/0.025
    rhoBetas = [rhoX_val*invT for rhoX_val in X]
    rho = Q[:]

    misccolors = ['#FF3232','#FB720E','#FC913A','#F9D423','#E6C049','#E1F5C4','#629E4F','#307B86','#21224A','#35353C', '#4D5E43', '#FFC132', '#FF7A32']
    misccolors = ['#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2','#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5']

    alphas = [0.001,0.1,1.0,5.0,8.0,10.0,20.0]
    betas  = list(np.linspace(0,15,201))

    #alphas = [0.001,1.0]
    #betas  = list(np.linspace(0,15,101))

    #betas  = list(np.linspace(0,15,100))
    #alphas = [0.001,0.1,1.0,5.0,8.0,10.0,20.0]

    #time1   = np.linspace(0,94,NT)
    #time1 = [float(i)*0.0188 for i in range(5001)]
    NT1 = 1e3
    NT2 = 5e3
    #NT2 = 1e5
    time1 = np.linspace(0,94,NT1)
    time2 = np.linspace(0,94,NT2)
    sab1,H,F = simpleGASKET(rhoBetas,rho,time1,alphas,betas)
    sab2,H,F = simpleGASKET(rhoBetas,rho,time2,alphas,betas)

    scalarMap = prepPlot2(alphas)

    for a in range(len(alphas)):
        plt.plot(betas[:-2],[sab1[b+a*len(betas)] for b in range(len(betas)-2)],\
                 color=scalarMap.to_rgba(a),linewidth=1.5,alpha=0.8)
        plt.plot(betas[:-2],[sab2[b+a*len(betas)] for b in range(len(betas)-2)],\
                 color=scalarMap.to_rgba(a),linewidth=1.5,alpha=0.8,linestyle='dashed')

    plt.show()


