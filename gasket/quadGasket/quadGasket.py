import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad
from scipy.interpolate import interp1d
from numpy import inf, cos, exp, tanh, sin


omega = ([0.0001,1.275E-2,1.9125E-2,2.55E-2,3.1875E-2,3.825E-2,4.4625E-2,5.1E-2,5.7375E-2, 6.375E-2, 6.63E-2, 6.885E-2, 7.14E-2, 7.395E-2, 7.65E-2, 8.2875E-2, 8.925E-2, 9.5625E-2, 1.02E-1, 1.08375E-1, 1.1475E-1, 1.21125E-1, 1.275E-1, 1.33875E-1, 0.15])

q = ([0.0, 2.25E-1,2.0E-1,8.125E-2,7.0E-2,7.125E-2,7.5E-2,7.9E-2,8.5E-2,9.5E-2, 1.15E-1, 1.197E-1, 1.214E-1, 1.218E-1, 1.195E-1, 1.125E-1, 9.75E-2, 8.71E-2, 7.91E-2, 7.35E-2, 6.88E-2, 6.5E-2, 6.1E-2, 5.71E-2,0.0])


def normalize(q,omega):
    integral = np.trapz(q,omega)
    am = 1.0/1.0086654 #convert to neutron mass
    norm1 = 1/(am*integral)
    q = [x*norm1 for x in q]

def get_H_F(omega,q,kbT,t):
    fterm = np.zeros(len(omega))
    hterm = np.zeros(len(omega))
    for i in range(len(omega)):
        fterm[i] = q[i]/omega[i]
        hterm[i] = fterm[i]/tanh(omega[i]/(2*kbT))
    
    
    Hgasket = np.zeros(len(t))
    Fgasket = np.zeros(len(t))
    for i in range(len(t)):
        sm = omega[0]*t[i]
        sinsm = math.sin(sm)
        cossm = math.cos(sm)
        for j in range(1,len(omega)):
            s = omega[j]*t[i]
            u = s-sm
            sins = math.sin(s)
            coss = math.cos(s)
            if u<0.005:
                st = u**2/6.0 - u**4/120.0
                ct = 0.5*u    - u**3/24.0
            else:
                sint = sins*cossm-coss*sinsm
                cost = coss*cossm+sins*sinsm
                st =  1.0-sint /u
                ct = (1.0-cost)/u
            Hgasket[i] = Hgasket[i] + hterm[j]*(st*sins+ct*coss)-hterm[j-1]*(st*sinsm-ct*cossm)
            Fgasket[i] = Fgasket[i] + fterm[j]*(ct*sins-st*coss)+fterm[j-1]*(ct*sinsm+st*cossm)
            sm    = s
            sinsm = sins
            cossm = coss
        Hgasket[i] = Hgasket[i]/t[i]
        Fgasket[i] = Fgasket[i]/t[i]
    
    return Fgasket,Hgasket





# tmin = 0.0001
# tmax = 15000.0
# alpha = 0.1
# beta = -13.5
# u = beta*am*kbT
# dbw = exp(-psq*Hgasket[0])
# Q1int,err1 = quad(qtable, tmin, tmax, weight='cos', wvar=beta*kbT, limit=500);
# Q2int,err1 = quad(q2table,tmin, tmax, weight='cos', wvar=beta*kbT, limit=500);
# Rint ,err2 = quad(rtable, tmin, tmax, weight='sin', wvar=beta*kbT, limit=500);



def quadGASKET(alphas,betas,kbT,omega,q,t):
    normalize(q,omega)
    Fgasket,Hgasket = get_H_F(omega,q,kbT,t)
    lambda_s = Hgasket[0]
    #plt.plot(t,Hgasket); plt.plot(t,Fgasket)
    #plt.show(); exit()
    ftable = interp1d(t,Fgasket,kind=5,fill_value=(0.0,0.0), bounds_error=False)
    htable = interp1d(t,Hgasket,kind=5,fill_value=(0.0,0.0), bounds_error=False)

    tmin = 0.0001; tmax = inf
    am = 1.0/1.0086654 #convert to neutron mass
    sab = [0.0]*len(alphas)*len(betas)

    for a,alpha in enumerate(alphas):
        print(alpha)
        psq = alpha*am*kbT
        #Sab = np.zeros(len(betas))
        psq = alpha*am*kbT
        dbw = exp(-psq*lambda_s)
        def q2table(t):
            return -dbw
        def qtable(t):
            return cos(psq*ftable(t))*exp(-psq*(lambda_s-htable(t)))
        def rtable(t):
            return sin(psq*ftable(t))*exp(-psq*(lambda_s-htable(t)))

        for b,beta in enumerate(betas):
            print(beta)
            Qint, err1 = quad(qtable, tmin, tmax, weight='cos', wvar=beta*kbT, epsrel=1e-1)#, limit=1500)
            Q2int,err1 = quad(q2table,tmin, tmax, weight='cos', wvar=beta*kbT, epsrel=1e-1)#, limit=1500)
            Rint, err2 = quad(rtable, tmin, tmax, weight='sin', wvar=beta*kbT, epsrel=1e-1)#, limit=1500)
            #Sab[b] = (Qint+Q2int-Rint)/math.pi
            sab[b+a*len(betas)] = (Qint+Q2int-Rint)/math.pi*np.exp(beta*0.5)
    return sab,Hgasket,Fgasket



if __name__=='__main__':
    alphas = [1e-5,1e-4,1e-3,1e-2,1e-1]
    betas = np.linspace(0,20.0,42)
    t = np.linspace(0.0001,15000.0, 5000)
    kbT = 0.0255 
    sab,H,F = quadGASKET(alphas,betas,kbT,omega,q,t)

"""

plt.semilogy(betas,Sab);

# T = 0.0255 #temperature in kT
temp = np.array([0.0255, 0.05, 0.075, 0.1, 0.125])
hterm_temp = np.zeros([len(omega),len(temp)])
for k in range(len(temp)):
    for i in range(len(omega)):
        z = omega[i]*0.5/temp[k]
        hterm_temp[i,k] = fterm[i]/tanh(z)




# How does Hgasket vary with T?

Hgasket_temp = np.zeros([len(t),len(temp)])
for k in range(len(temp)):
    for i in range(len(t)):
        sm = omega[0]*t[i]
        sinsm = math.sin(sm)
        cossm = math.cos(sm)
        for j in range(1,len(omega)):
            s = omega[j]*t[i]
            u = s-sm
            sins = math.sin(s)
            coss = math.cos(s)
            if u<0.005:
                st = u**2/6.0 - u**4/120.0
                ct = 0.5*u - u**3/24.0
            else:
                sint = sins*cossm-coss*sinsm
                cost = coss*cossm+sins*sinsm
                st = 1.0-sint/u
                ct = (1.0-cost)/u
            Hgasket_temp[i,k] = Hgasket_temp[i,k] + hterm_temp[j,k]*(st*sins+ct*coss)-hterm_temp[j-1,k]*(st*sinsm-ct*cossm)
            sm = s
            sinsm=sins
            cossm=coss
        Hgasket_temp[i,k] = Hgasket_temp[i,k]/t[i]



plt.plot(t[0:1000],Hgasket_temp[0:1000,0])
plt.plot(t[0:1000],Hgasket_temp[0:1000,1])
plt.plot(t[0:1000],Hgasket_temp[0:1000,2])
plt.plot(t[0:1000],Hgasket_temp[0:1000,3])
plt.plot(t[0:1000],Hgasket_temp[0:1000,4])
"""

