import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import cos, sin, exp
from math import pi


#0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0
# Define and normalize that phonon distribution (two oscillators wit combined 
# weight of 0.5, with contin weight of 0.5)
#0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0
x = ([6.375E-3, 1.275E-2, 1.9125E-2, 2.55E-2, 3.1875E-2, 3.825E-2, 4.4625E-2, 
      5.1E-2, 5.7375E-2, 6.375E-2, 6.63E-2, 6.885E-2, 7.14E-2, 7.395E-2, 
      7.65E-2, 8.2875E-2, 8.925E-2, 9.5625E-2, 1.02E-1, 1.08375E-1, 1.1475E-1, 
      1.21125E-1, 1.275E-1, 1.33875E-1, 1.4025E-1, 1.46625E-1, 1.53E-1, 
      1.59375E-1, 1.6575E-1, 1.75E-1, 0.204, 0.205, 0.206, 0.479, 0.480, 0.481])

rho = ([2.25E-1, 2.0E-1, 8.125E-2, 7.0E-2, 7.125E-2, 7.5E-2, 7.9E-2, 8.5E-2, 
      9.5E-2, 1.15E-1, 1.197E-1, 1.214E-1, 1.218E-1, 1.195E-1, 1.125E-1, 9.75E-2, 
      8.71E-2, 7.91E-2, 7.35E-2, 6.88E-2, 6.5E-2, 6.1E-2, 5.71E-2, 5.4E-2, 
      5.15E-2, 4.88E-2, 4.59E-2, 4.31E-2, 4.2E-2, 1.0E-10, 1.0E-10, 0.333, 
      1.0E-10, 1.0E-10, 0.667, 1.0E-10 ])

integral = np.trapz(rho[0:28],x[0:28]) + rho[0]*x[0]/3
am = 1.0/1.0086654 #convert to neutron mass
continWgt = 0.5
discreWgt = 0.5
norm1 = continWgt/(am*integral)
norm2 = discreWgt/am/((x[31]-x[30])*rho[31]+(x[34]-x[33])*rho[34])
rho = [x*norm1 for x in rho[:29]] + [x*norm2 for x in rho[29:]]

#0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0





total_time = 5000
time_steps = 5000
time = np.linspace(0.0,total_time,time_steps)


T = 0.0255 #temperature in kT
integrand, coth = np.zeros(len(rho)), np.zeros(len(rho))
for i in range(len(rho)):
    coth[i] = math.cosh(x[i]*0.5/T)/math.sinh(x[i]*0.5/T)
    integrand[i] = rho[i]/x[i]*coth[i]


integral = np.trapz(integrand,x)
invF  = 2.0*T/x[0]

F = np.zeros(time_steps)
H = np.zeros(time_steps)
H[0] = 0.5*rho[0]*(invF+coth[0]) + integral








#0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0
# Prepare the F and H functions. 
#0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0
for t in range(1,time_steps):
    u = x[0]*time[t]

    c  = (1-cos(u))/u           if u >= 0.005 else 0.5*u+u**3/24.0
    s  =    sin(u)              if u >= 0.005 else     u-u**3/6.0
    cs = sin(u)/u**2 - cos(u)/u if u >= 0.005 else u/3.0 - u**3/30.0

    F[t] += time[t]*rho[0]*cs
    H[t] += time[t]*rho[0]/u*(coth[0]*(s-c)+c*invF)

    sinLast = sin(x[0]*time[t])
    cosLast = cos(x[0]*time[t])

    for i in range(1,len(rho)):
        sinNext = sin(x[i]*time[t])
        cosNext = cos(x[i]*time[t])

        theta = x[i]*time[t]-x[i-1]*time[t]
        st = 1-sin(theta) /theta if theta >= 5e-3 else theta**2/6 - theta**4/120
        ct =(1-cos(theta))/theta if theta >= 5e-3 else 0.5*theta - theta**3/24

        H[t] += rho[i]/x[i]    *coth[i]  *(st*sinNext + ct*cosNext) - \
                rho[i-1]/x[i-1]*coth[i-1]*(st*sinLast  - ct*cosLast )

        F[t] += rho[i]  /x[i]  *(ct*sinNext-st*cosNext) + \
                rho[i-1]/x[i-1]*(ct*sinLast+st*cosLast)

        sinLast = sinNext
        cosLast = cosNext

    H[t] /= time[t]
    F[t] /= time[t]

#0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0Oo.oO0








def gasket(alpha,beta):
    betaT = beta*T #AL in fortran code
    psq = alpha*am*T  # B in fortran code
    dbw = math.exp(-psq*H[0])  
    #dbw2 = math.exp(-psq*H2[0])  
    #dbwp = dbw/math.pi   #F in fortran code
    S = 0.0
    sm = time[0]*betaT
    sinLast = sin(sm)
    cosLast = cos(sm)
    #ex1 = math.exp(psq*H[0])
    #qLast = (math.cos(psq*F[0])*ex1 -1.0)
    #r1 = math.sin(psq*F[0])*ex1
    ex1 = 1.0
    qLast = cos(psq*F[0])*ex1 - dbw
    rLast = sin(psq*F[0])*ex1
    for t in range(1,time_steps):
        sinNext = sin(time[t]*betaT)
        cosNext = cos(time[t]*betaT)

        theta = time[t]*betaT-time[t-1]*betaT

        st = 1-sin(theta) /theta if theta >= 5e-3 else theta**2/6 - theta**4/120
        ct =(1-cos(theta))/theta if theta >= 5e-3 else 0.5*theta - theta**3/24

        #ex1 = math.exp(psq*H[t])
        #qNext = (math.cos(psq*F[t])*ex1 - 1.0)
        #rNext = math.sin(psq*F[t])*ex1
        ex1   = exp(-psq*(H[0]-H[t]))
        qNext = cos(psq*F[t])*ex1 - dbw
        rNext = sin(psq*F[t])*ex1

        S += qNext*(st*sinNext+ct*cosNext) - qLast*(st*sinLast-ct*cosLast) - \
             rNext*(ct*sinNext-st*cosNext) - rLast*(st*cosLast+ct*sinLast)
        sinLast = sinNext
        cosLast = cosNext
        qLast = qNext
        rLast = rNext
    #S=S*dbwp/betaT
    return S/(betaT*pi)



beta = 0.5
Sasy = gasket(0.5,0.5)
Ssym = Sasy/math.exp(beta/2.)

x = gasket(5e-3,5e-3)
if abs(x-0.6616955610028908)/x > 1e-6:
    print(1,x,0.6616955610028908)
    exit()

x = gasket(0.5,0.5)
if abs(x-8.286163606448213)/x> 1e-6:
    print(2,x,8.286163606448213)
    exit()

x = gasket(1.5,1.5)
if abs(x-1.042763620113287 )/x  > 1e-6: 
    print(3,x,1.042763620113287)
    exit()

x = gasket(3.0,3.0) 
if abs(x-0.1616246706325477)/x  > 1e-6:
    print(4,x,0.1616246706325477)
    exit()

x = gasket(1.0,8.0)
if abs(x-0.00022785000823275238)/x > 1e-6:
    print(5,x,0.00022785000823275238)
    exit()

x = gasket(8.0,1.0)
if abs(x-1.2160162511496266)/x > 1e-6:
    print(6,x,1.2160162511496266)
    exit()

#x = gasket(15.0,15.0)
#if abs(x-4.677122999348168e-07)/x > 1e-6:
#    print(7,x,4.677122999348168e-07)
#    exit()

print("All good!")


"""

# In[13]:


sigmab = 10.0
A = 1.0
E = 0.1
kT = 0.025
kTeff = 0.035
Emin = max(0.0, E - 16*kT)
Emax = E + 16*kT
Eprime = np.linspace(Emin,Emax,200)
mu = np.linspace(-1.0,1.0,10000)


# In[14]:


sigma = np.zeros([len(Eprime),len(mu)])
sigma_sct = np.zeros([len(Eprime),len(mu)])
sigma_gt = np.zeros([len(Eprime),len(mu)])
for i in range(len(Eprime)):
    #print(i)
    for j in range(len(mu)):
        alpha = (Eprime[i]+E - 2*mu[j]*math.sqrt(Eprime[i]*E))/A/kT
        beta = (Eprime[i]-E)/kT
        if alpha<0.001:
            alpha=0.001
        #free gas
        S = 1/math.sqrt(4*math.pi*alpha)*math.exp(-((alpha+beta)**2)/4/alpha)
        #sct
        St = 1/math.sqrt(4*math.pi*alpha*kTeff/kT)*math.exp(-((alpha+beta)**2)/4/(alpha*kTeff/kT))
        #Gasket
        Sg = gasket(alpha,-beta)
        sigma[i,j] = sigmab/2/kT*math.sqrt(Eprime[i]/E)*S
        sigma_sct[i,j] = sigmab/2/kT*math.sqrt(Eprime[i]/E)*St
        sigma_gt[i,j] = sigmab/2/kT*math.sqrt(Eprime[i]/E)*Sg


# In[15]:


#angle integrated E' prime
sigmaE = np.sum(sigma, axis=1)
sigmaE_sct = np.sum(sigma_sct, axis=1)
sigmaE_gt = np.sum(sigma_gt, axis=1)
sigmaE_norm = sigmaE/np.linalg.norm(sigmaE)
sigmaE_sct_norm = sigmaE_sct/np.linalg.norm(sigmaE_sct)
sigmaE_gt_norm = sigmaE_gt/np.linalg.norm(sigmaE_gt)


# In[16]:


plt.plot(Eprime, sigmaE_norm,'b')
plt.plot(Eprime, sigmaE_sct_norm,'r')
plt.plot(Eprime, sigmaE_gt_norm,'g')


# In[46]:


E = 1.0
kT = 0.025
kTeff = 0.035
A = 1.0
betamin = E/kT
beta = np.linspace(-betamin,20.0,45)
#alpha = np.linspace(0.1,200,200)
S = np.zeros(len(beta))
St = np.zeros(len(beta))
Sg = np.zeros(len(beta))
for i in range(len(beta)):
    #print(i)
    alphamin = (math.sqrt(E)-math.sqrt(E+beta[i]*kT))**2/A/kT
    alphamax = (math.sqrt(E)+math.sqrt(E+beta[i]*kT))**2/A/kT
    if alphamin < 0.1:
        alphamin = 0.1
    #if alphamax > 200.0:
    #    alphamax = 200.0
    #print(alphamin, alphamax)
    alpha = np.linspace(alphamin,alphamax,200)
    for j in range(len(alpha)):
        #free gas
        #S[i] = S[i] + 1/math.sqrt(4*math.pi*alpha[j])*math.exp(-((alpha[j]+beta[i])**2)/4/alpha[j])
        #sct
        #St[i] = St[i] + 1/math.sqrt(4*math.pi*alpha[j]*kTeff/kT)*math.exp(-((alpha[j]+beta[i])**2)/4/(alpha[j]*kTeff/kT))
        #Gasket
        if alpha[j] < 30.0:
            Sg[i] = Sg[i] + gasket(alpha[j],beta[i])
        else:
            Sg[i] = Sg[i] + 1/math.sqrt(4*math.pi*alpha[j]*kTeff/kT)*math.exp(-((alpha[j]+beta[i])**2)/4/(alpha[j]*kTeff/kT))
        


# In[47]:


#plt.plot(beta,S)
#plt.plot(beta, St)
plt.plot(beta,Sg)


# In[28]:


Sg

"""
