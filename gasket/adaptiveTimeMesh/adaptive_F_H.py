import sys
sys.path.append('../../phononDistributions')
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,exp,sinh,cosh
from waterData import *


# What this is meant to show is that converging your t grid to represent F(t)
# and H(t) accurately is not realy a worthwhile endeavor, since it actually 
# makes the final S(a,b) answer way worse (while using more points) than the
# typical uniform approach.


def get_CT(x):  return (1.0-cos(x))/x       if x > 5e-3 else x*0.5 - x**3/24.0
def get_ST(x):  return  1.0-sin(x) /x       if x > 5e-3 else x**2/6. - x**4/120.
def get_sin(x): return sin(x)               if x > 5e-3 else x - x**3*0.1666666
def get_F0(x):  return sin(x)/x**2-cos(x)/x if x > 5e-3 else x/3. - (x**3)/30.

rho      = Q[:]
betas = [x/0.0255 for x in X]
coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in betas]

def getFtHt(t,rho,betas,norm):
    assert(t>0.0)

    C = get_CT( betas[0]*t)
    S = get_sin(betas[0]*t)
    Ht = rho[0]/(betas[0])*(coth[0]*(S-C)+C/(betas[0]*0.5))
    Ft = rho[0]*t*get_F0(betas[0]*t)
    sinNext = sin(betas[0]*t)
    cosNext = cos(betas[0]*t)

    for j in range(1,len(betas)):
        sinLast = sinNext; sinNext = sin(betas[j]*t) # Advance sin
        cosLast = cosNext; cosNext = cos(betas[j]*t) # and cos
        if abs(betas[j]/betas[j-1] - 1.0) > 5e-7:
            theta = (betas[j]-betas[j-1])*t;
            ST, CT = get_ST(theta), get_CT(theta)
        Ht += rho[j  ]/betas[j  ] * (ST*sinNext + CT*cosNext ) * coth[j  ] - \
              rho[j-1]/betas[j-1] * (ST*sinLast - CT*cosLast ) * coth[j-1]
        Ft += rho[j  ]/betas[j  ] * (CT*sinNext - ST*cosNext ) + \
              rho[j-1]/betas[j-1] * (CT*sinLast + ST*cosLast )
    return Ft*norm/t,Ht*norm/t



def get_F_H_uniform(betas,Q,time,norm):
    coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in betas]
    F, H = [0.0]*len(time), [0.0]*len(time)
    norm = 1.0 / (Q[0]*betas[0]/3.0 + np.trapz(Q,betas))

    # Putting in the t=0 terms
    F[0] = 0.0;
    H[0] = norm * ( Q[0]/betas[0] + 0.5*Q[0]*coth[0] +     \
                    np.trapz( [Q[i]*coth[i]/betas[i] for i in range(len(Q))], betas ) )

    for t in range(1,len(time)):
        F[t],H[t] = getFtHt(time[t],Q,betas,norm)

    return F,H


def get_F_H_adaptively(rho,betas,spacing,max_t_val,tolerance,norm):
    xL = 0.0; final_X = [xL]
    # F0 = 0, H0 = something unfortunate (these are the t=0 contributions)
    yL = (0.0, norm*(rho[0]/betas[0] + 0.5*rho[0]*coth[0] + \
               np.trapz([rho[i]*coth[i]/betas[i] \
               for i in range(len(betas))],betas)))

    xR = xL+0.1; yR = getFtHt(xR,rho,betas,norm)
    final_F = [yL[0]]
    final_H = [yL[1]]
    # Set convergence criteria. So I'm going to look at the last N values and 
    # make sure they're all less than ratio*max(F) or ratio*max(H)
    N = 10; ratio = 1e-5
    max_vals = (0.0,0.0)

    while xR <= max_t_val:
        xM = (xL+xR)*0.5
        yM = getFtHt(xM,rho,betas,norm)
        fM,hM = yM
        guessF = (yL[0]+yR[0])*0.5
        guessH = (yL[1]+yR[1])*0.5
        error1 = abs(guessF-yM[0]) if yM[0] < 1e-5 else abs((guessF-yM[0])/yM[0])
        error2 = abs(guessH-yM[1]) if yM[1] < 1e-5 else abs((guessH-yM[1])/yM[1]) 
        if error1 < tolerance and (xL == 0.0 or error2 < tolerance):
            final_X.append(xR)
            final_F.append(yR[0])
            final_H.append(yR[1])
            max_vals = ( max(abs(yR[0]),max_vals[0]), \
                         max(abs(yR[1]),max_vals[1]) )
            xL = xR; xR = xL + spacing
            yL = yR; yR = getFtHt(xR,rho,betas,norm)

            if len(final_X) > N:
                if max([abs(x) for x in final_F[-N:]]) < ratio*max_vals[0] and \
                   max([abs(x) for x in final_H[-N:]]) < ratio*max_vals[1]:
                    return final_X,final_F,final_H
        else:
            xR = xM; yR = yM
    print("Didn't converge, maybe try increasing max t value")
    return final_X,final_F,final_H

def plot_sab(time,F,H,a_s,b_s,colors,linestyle):
    alpha_exp = [exp(-alpha*H[0]) for alpha in a_s]
    beta_exp  = [exp(-beta*0.5)   for beta  in b_s]
    for a,alpha in enumerate(a_s):
        alpha_H_exp = [exp(alpha*H[t]) for t in range(len(time))]
        sin_alpha_F = [sin(alpha*F[t]) for t in range(len(time))]
        cos_alpha_F = [cos(alpha*F[t]) for t in range(len(time))]
        sabChunk = [0.0]*len(b_s)
        for b,beta in enumerate(b_s):
            R = [alpha_H_exp[t] * cos_alpha_F[t] * cos(beta*time[t]) \
               - alpha_H_exp[t] * sin_alpha_F[t] * sin(beta*time[t]) \
               - cos(beta*time[t]) for t in range(len(time))]
            sabChunk[b] = np.trapz(R,time) * alpha_exp[a] * beta_exp[b]
        plt.plot(b_s,sabChunk,label=str(alpha),color=colors[a],linestyle=linestyle,linewidth=2)





rhoBetas = [x/0.0255 for x in X]
rho = Q[:]
invArea = 1.0/np.trapz(rho,rhoBetas); rho = [invArea * x for x in rho]
norm = 1.0 / (rho[0]*rhoBetas[0]/3.0 + np.trapz(rho,rhoBetas))


a_s = [0.1,1.0,5.0,10.0,20.0]
b_s = np.linspace(0,20,100)
colors = ["#F8B195", "#F67280", "#C06C84", "#6C5B7B", "#355C7D"]
colors = ['#581845','#900C3F','#C70039','#FF5733','#FFC30F']

spacing = 0.1; max_t_val = 300.0; tolerance = 1e-3
time = np.linspace(0,94,941)
F_normal,   H_normal   = get_F_H_uniform(rhoBetas,rho,time,norm)
time_adaptive, F_adaptive, H_adaptive = get_F_H_adaptively(rho,rhoBetas,\
     spacing=time[1]-time[0],max_t_val=max_t_val,tolerance=1e-3,norm=norm)

print("Constructed F(t) and H(t) using the normal, uniform grid")
print("with",len(time),'values where tMax = ',time[-1],'\n')
print("Constructed F(t) and H(t) using an adaptive time grid")
print("with",len(time_adaptive),'values where tMax = ',time_adaptive[-1])

wantToPlotFH  = True
wantToPlotSab = True

if wantToPlotFH:
    fig = plt.figure(num=None,figsize=(12,6),dpi=100,facecolor='w',edgecolor='k') 
    plt.rcParams.update({'font.size': 14})
    plt.plot(time,F_normal,label='F(t) (uniform)',color='r',linestyle='solid'); 
    plt.plot(time,H_normal,label='H(t) (uniform)',color='b',linestyle='solid'); 
    plt.plot(time_adaptive,F_adaptive,label='F(t) (adaptive)',color='r',linestyle='dashed'); 
    plt.plot(time_adaptive,H_adaptive,label='F(t) (adaptive)',color='b',linestyle='dashed'); 
    plt.title('F(t) and H(t) for Uniform and Adaptive Time Grids')
    plt.legend(loc='best'); plt.xlabel('time'); plt.show()

if wantToPlotSab:
    fig = plt.figure(num=None,figsize=(15,6),dpi=100,facecolor='w',edgecolor='k') 
    plt.rcParams.update({'font.size': 14})
    plt.subplot(121)
    plt.title('S(a,-b) for uniform t grid'); plt.ylabel('S(a,-b)'); plt.xlabel('beta')
    plot_sab(time         ,F_normal  ,H_normal  ,a_s,b_s,colors,'solid'); plt.yscale('log')
    plt.legend(loc='best'); 
    plt.subplot(122)
    plt.title('S(a,-b) for adaptive t grid'); plt.ylabel('S(a,-b)'); plt.xlabel('beta')
    plot_sab(time_adaptive,F_adaptive,H_adaptive,a_s,b_s,colors,'solid'); plt.yscale('log')
    plt.legend(loc='best'); 
    plt.show()


