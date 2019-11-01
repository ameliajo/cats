import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.colors as colors
import numpy as np 
from numpy import sin,cos,sinh,cosh,exp
from math import pi
from get_F_H import *
import math
from scipy.integrate import quad
from scipy.interpolate import interp1d
from numpy import inf, cos, exp, tanh, sin



def prepPlot(alphas):
    figure(num=None, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    mp.rcParams.update({'font.size': 15})
    cnorm = mp.colors.Normalize(vmin=0,vmax=len(alphas))
    scalarMap = mp.cm.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) 
    colorList = [mp.colors.rgb2hex(scalarMap.to_rgba(i)) for i in range(len(alphas))]
    colorBar = plt.contourf([[0,0],[0,0]], list(range(len(alphas)+1)), \
               cmap=mp.colors.ListedColormap(colorList))
    plt.clf()

    plt.title('H in H2O Scattering Law (Direct Method)')
    plt.xlabel('beta'); plt.ylabel('S(a,b) (Non-symmetric)'); plt.yscale('log')

    cbar = plt.colorbar(colorBar);     cbar.ax.get_yaxis().set_ticks([])
    cbar.ax.get_yaxis().labelpad = 40; cbar.ax.set_ylabel('alpha')
    for j, lab in enumerate(["    "+str(a) for a in alphas[:]]):
        cbar.ax.text(1.7,(2*j+1)/len(alphas), lab, ha='center', va='center')
    return scalarMap



def simpleGASKET(rhoBetas,rho,time,alphas,betas):
    print("Simple GASKET")
    print("-------------")

    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho     = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)

    tmin, tmax = 0, time[-1]
    ftable = interp1d(time,F,kind=0,fill_value=(0.0,0.0), bounds_error=False)
    htable = interp1d(time,H,kind=0,fill_value=(0.0,0.0), bounds_error=False)

    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0]) for alpha in alphas] # H[0] = debye waller 
    beta_exp  = [exp(beta*0.5)    for beta  in betas ] # This is to turn 
                                         # the scattering law to be symmetric
    for a,alpha in enumerate(alphas):
        print(alpha)
        alpha_H_exp = [exp(alpha*htable(time[t])) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*ftable(time[t])) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*ftable(time[t])) for t in range(len(time))]
        for b,beta in enumerate(betas):
            print('     ',beta)
            integrand = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
                       - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
                       - cos(beta*time[t]) for t in range(len(time))]
            sab[b+a*len(betas)] = np.trapz(integrand,time) * alpha_exp[a] * beta_exp[b] / math.pi
    return sab,H,F




def quadGASKET(rhoBetas,rho,time,alphas,betas):
    print("Quad GASKET")
    print("-----------")

    invArea = 1.0/np.trapz(rho,rhoBetas)
    rho     = [invArea * x for x in rho]
    F, H    = get_F_H(rhoBetas,rho,time)

    tmin, tmax = 0, time[-1]
    ftable = interp1d(time,F,kind=1,fill_value=(0.0,0.0), bounds_error=False)
    htable = interp1d(time,H,kind=1,fill_value=(0.0,0.0), bounds_error=False)

    sab       = [0.0]*len(betas)*len(alphas)
    alpha_exp = [exp(-alpha*H[0]) for alpha in alphas] # H[0] = debye waller 
    beta_exp  = [exp(beta*0.5)    for beta  in betas ] # This is to turn 
                                         # the scattering law to be symmetric
    lambda_s = H[0]
    getDbw = lambda t,alpha: -exp(-alpha*lambda_s)
    getQ   = lambda t,alpha:  cos(alpha*ftable(t))*exp(-alpha*(lambda_s-htable(t)))
    getR   = lambda t,alpha:  sin(alpha*ftable(t))*exp(-alpha*(lambda_s-htable(t)))

    for a,alpha in enumerate(alphas):
        print(alpha)
        alpha_H_exp = [exp(alpha*htable(time[t])) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*ftable(time[t])) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*ftable(time[t])) for t in range(len(time))]
        for b,beta in enumerate(betas):
            print('     ',beta)
            Qint, err1 = quad(getQ,  tmin, tmax, args=(alpha,), weight='cos', wvar=beta)#, epsabs=1e-5)
            Q2int,err1 = quad(getDbw,tmin, tmax, args=(alpha,), weight='cos', wvar=beta)#, epsabs=1e-5)
            Rint, err2 = quad(getR,  tmin, tmax, args=(alpha,), weight='sin', wvar=beta)#, epsabs=1e-5)
            sab[b+a*len(betas)] = (Qint+Q2int-Rint)/math.pi * beta_exp[b]

    return sab,H,F





if __name__=="__main__":
    import sys
    sys.path.append("../phononDistributions"); from waterDataContinuous import *

    invT = 1.0/0.025
    rhoBetas = [rhoX_val*invT for rhoX_val in X]

    misccolors = ['#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a',\
        '#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2',\
        '#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5']

    alphas = [1e-5,1e-4,1e-3,1e-2,0.1]#,1.0,5.0,8.0,10.0,20.0]
    betas  = list(np.linspace(0,10,41))

    NT1 = 1e3; time1 = np.linspace(0,94,NT1)

    sab1,H,F = simpleGASKET(rhoBetas,Q,time1,alphas,betas)
    sab2,H,F = quadGASKET(rhoBetas,Q,time1,alphas,betas)

    scalarMap = prepPlot(alphas);  plt.clf()
    colorList = [colors.rgb2hex(scalarMap.to_rgba(i)) for i in range(len(alphas))]
    mymap = colors.ListedColormap(colorList)
    colorBar = plt.contourf([[0,0],[0,0]], list(range(len(alphas)+1)), cmap=mymap)
    plt.clf()
    cbar = plt.colorbar(colorBar); cbar.ax.get_yaxis().set_ticks([])
    for j, lab in enumerate(["    "+str(a) for a in alphas[:]]):
        cbar.ax.text(1.7,(2*j+1)/(2*len(alphas)),lab,ha='center',va='center')

    for a in range(len(alphas)):
        plt.plot(betas[:-2],[sab1[b+a*len(betas)] for b in range(len(betas)-2)],\
                 color=scalarMap.to_rgba(a),linewidth=1.5,alpha=0.8)
        plt.plot(betas[:-2],[sab2[b+a*len(betas)] for b in range(len(betas)-2)],\
                 color=scalarMap.to_rgba(a),linewidth=1.5,alpha=0.8,linestyle='dashed')

    plt.yscale('log')
    plt.show()


