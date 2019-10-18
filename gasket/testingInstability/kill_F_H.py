import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,exp,sinh,cosh
import sys
sys.path.append('../phononDistributions')
from waterDataContinuous import *

def get_CT(x): return (1.-cos(x))/x          if x > 5e-3 else x*0.5 - x**3/24.
def get_ST(x): return  1.-sin(x) /x          if x > 5e-3 else x**2/6. - x**4/120.
def get_sin(x):return sin(x)                 if x > 5e-3 else x - x**3*0.1666666
def get_F0(x): return sin(x)/x**2 - cos(x)/x if x > 5e-3 else x/3. - (x**3)/30.

energies = X[:]
rho      = Q[:]
betas = [x/0.0255 for x in energies]
coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in betas]


def get_Ft_Ht(betas,Q,t,norm):
    C = get_CT(betas[0]*time[t])
    S = get_sin(betas[0]*time[t])
    Ht = Q[0]/(betas[0])*(coth[0]*(S-C)+C/(betas[0]*0.5));
    Ft = Q[0]*time[t]*get_F0(betas[0]*time[t])

    sinNext = sin(betas[0]*time[t])
    cosNext = cos(betas[0]*time[t])

    for j in range(1,len(betas)):
        sinLast = sinNext; sinNext = sin(betas[j]*time[t]) # Advance sin
        cosLast = cosNext; cosNext = cos(betas[j]*time[t]) # and cos

        if abs(betas[j]/betas[j-1] - 1.0) > 5e-7:
            theta = (betas[j]-betas[j-1])*time[t];
            ST, CT = get_ST(theta), get_CT(theta)

        Ht += Q[j  ]/betas[j  ] * (ST*sinNext + CT*cosNext ) * coth[j  ] - \
              Q[j-1]/betas[j-1] * (ST*sinLast - CT*cosLast ) * coth[j-1]
        Ft += Q[j  ]/betas[j  ] * (CT*sinNext - ST*cosNext ) + \
              Q[j-1]/betas[j-1] * (CT*sinLast + ST*cosLast )

    return Ft*norm/time[t], Ht*norm/time[t]


def get_F_H(betas,Q,time,norm):
    coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in betas]
    F, H = [0.0]*len(time), [0.0]*len(time)
    norm = 1.0 / (Q[0]*betas[0]/3.0 + np.trapz(Q,betas))

    # Putting in the t=0 terms
    F[0] = 0.0;
    H[0] = norm * ( Q[0]/betas[0] + 0.5*Q[0]*coth[0] +     \
                    np.trapz( [Q[i]*coth[i]/betas[i] for i in range(len(Q))], betas ) )

    for t in range(1,len(time)):
        F[t],H[t] = get_Ft_Ht(betas,Q,t,norm)
    return F,H




invT = 1.0/0.0255
rhoBetas = [x*invT for x in X]
rho = Q[:]
invArea = 1.0/np.trapz(rho,rhoBetas)
rho = [invArea * x for x in rho]
norm = 1.0 / (rho[0]*rhoBetas[0]/3.0 + np.trapz(rho,rhoBetas))

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


fig = plt.figure(num=None,figsize=(10,10),dpi=100,facecolor='w',edgecolor='k') 
a_s = [0.1,1.0,5.0,10.0,20.0]
b_s = np.linspace(0,20,100)
colors = ["#F8B195", "#F67280", "#C06C84", "#6C5B7B", "#355C7D"]
colors = ['#581845','#900C3F','#C70039','#FF5733','#FFC30F']

time = np.linspace(0,100,2001)
print('Youre trying this with a randomly chosen time grid - see')
print('what happens if you use a more carefully chosen time grid')
#time = np.linspace(0,96.14,941)
F, H    = get_F_H(rhoBetas,rho,time,norm)


plot_sab(time,F,H,a_s,b_s,colors,'solid'); plt.yscale('log')
F[-1] = 0.0; F.append(0.0)
H[-1] = 0.0; H.append(0.0)
time = list(time); time.append(time[-1]+(time[-1]-time[-2]))
plot_sab(time,F,H,a_s,b_s,colors,'dashed'); plt.yscale('log')
plt.title('This shows that adding a few zeros at end doesnt help')
plt.show()


tMid = int(len(time)*0.5)
scaling = [1.0]*tMid + \
          [exp(-(time[i]-time[tMid])**2/(0.5*(time[-1]-time[tMid]))) \
           for i in range(tMid,len(time))]

F2 = [F[i]*scaling[i] for i in range(len(F))]
H2 = [H[i]*scaling[i] for i in range(len(H))]

tMid = int(len(time)*0.8)
scaling = [1.0]*tMid + \
          [exp(-(time[i]-time[tMid])**2/(0.5*(time[-1]-time[tMid]))) \
           for i in range(tMid,len(time))]

F3 = [F[i]*scaling[i] for i in range(len(F))]
H3 = [H[i]*scaling[i] for i in range(len(H))]

plt.plot(time,F,label='F(t)',color='r',linestyle='solid'); 
plt.plot(time,H,label='H(t)',color='b',linestyle='solid'); 
plt.plot(time,F2,label='F2(t)',color='r',linestyle='dashed'); 
plt.plot(time,H2,label='H2(t)',color='b',linestyle='dashed'); 
plt.plot(time,F3,label='F3(t)',color='r',linestyle='dotted'); 
plt.plot(time,H3,label='H3(t)',color='b',linestyle='dotted'); 
plt.title('Zoom into the tail end - Im damping the F(t) and H(t)')
plt.show()


plot_sab(time,F,H,a_s,b_s,colors,'solid'); plt.yscale('log')
plot_sab(time,F2,H2,a_s,b_s,colors,'dashed'); plt.yscale('log')
plot_sab(time,F3,H3,a_s,b_s,colors,'dotted'); plt.yscale('log')
plt.title('Look at the effect this damping has on S(a,-b)')
plt.show()

