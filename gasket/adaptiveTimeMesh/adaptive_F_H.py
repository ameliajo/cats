import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,exp,sinh,cosh
from waterDataContinuous import *

def get_CT(x):
    return (1.0-cos(x))/x if x > 0.005 else x*0.5 - x**3/24.0
def get_ST(x):
    return  1.0-sin(x) /x if x > 0.005 else x**2/6.0 - x**4/120.0
def get_sin(x):
    return sin(x) if x > 0.005 else x - x**3*0.1666666
def get_F0(x):
    return sin(x)/x**2 - cos(x)/x if x > 0.005 else x/3.0 - (x**3)/30.0

energies = X[:]
rho      = Q[:]
betas = [x/0.0255 for x in energies]
coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in betas]

norm = 1.0 / (Q[0]*betas[0]/3.0 + np.trapz(Q,betas))

def getFt(t,rho,betas):
    assert(t>0.0)

    C = get_CT( betas[0]*t)
    S = get_sin(betas[0]*t)
    Ht = Q[0]/(betas[0])*(coth[0]*(S-C)+C/(betas[0]*0.5))
    Ft = rho[0]*t*get_F0(betas[0]*t)
    sinNext = sin(betas[0]*t)
    cosNext = cos(betas[0]*t)

    for j in range(1,len(betas)):
        sinLast = sinNext; sinNext = sin(betas[j]*t) # Advance sin
        cosLast = cosNext; cosNext = cos(betas[j]*t) # and cos
        if abs(betas[j]/betas[j-1] - 1.0) > 5e-7:
            theta = (betas[j]-betas[j-1])*t;
            ST, CT = get_ST(theta), get_CT(theta)
        Ht += Q[j  ]/betas[j  ] * (ST*sinNext + CT*cosNext ) * coth[j  ] - \
              Q[j-1]/betas[j-1] * (ST*sinLast - CT*cosLast ) * coth[j-1]
        Ft += rho[j  ]/betas[j  ] * (CT*sinNext - ST*cosNext ) + \
              rho[j-1]/betas[j-1] * (CT*sinLast + ST*cosLast )
    return Ft*norm/t,Ht*norm/t

def getF(rho,betas,spacing,max_t_val,tolerance):
    xL = 0.0; final_X = [xL]

    # Figure out what t = 0 contributions will be
    F0 = 0.0
    H0 = Q[0]/betas[0] + 0.5*Q[0]*coth[0] + np.trapz([Q[i]*coth[i]/betas[i] for i in range(len(betas))],betas)
    yL = (F0,H0)

    xR = xL+1e-7; yR = getFt(xR,rho,betas)
    final_F = [yL[0]]
    final_H = [yL[1]]
    max_F = 0.0
    max_H = 0.0
    while xR < max_t_val:
        xM = (xL+xR)*0.5
        yM = getFt(xM,rho,betas)
        guessF = (yL[0]+yR[0])*0.5
        guessH = (yL[1]+yR[1])*0.5
        error1 = abs(guessF-yM[0]) if yM[0] < 1e-5 else abs((guessF-yM[0])/yM[0])
        error2 = abs(guessH-yM[1]) if yM[1] < 1e-5 else abs((guessH-yM[1])/yM[1]) 
        if error1 < tolerance and (xL == 0.0 or error2 < tolerance):
            final_X.append(xR)
            final_F.append(yR[0])
            final_H.append(yR[1])
            if abs(yR[0]) > max_F: max_F = abs(yR[0])
            if abs(yR[1]) > max_H: max_H = abs(yR[1])
            xL = xR; xR = xL + spacing
            yL = yR; yR = getFt(xR,rho,betas)
            if len(final_X) > 20:
                if max([abs(x) for x in final_F[-20:]]) < 1e-3*max_F and \
                   max([abs(x) for x in final_H[-20:]]) < 1e-3*max_H:
                    return final_X,final_F,final_H
        else:
            xR = xM; yR = yM
    return final_X,final_F,final_H

spacing = 2.0
max_t_val = 100.0
tolerance = 1e-3
final_X, final_F, final_H = getF(rho,betas,spacing,max_t_val,tolerance)

inv_area = 1.0/np.trapz(final_F,final_X); final_F = [inv_area * F for F in final_F]
inv_area = 1.0/np.trapz(final_H,final_X); final_H = [inv_area * H for H in final_H]


plt.plot(final_X,final_F)
plt.plot(final_X,final_H)
#plt.plot(final_X,contrib_F)
#plt.plot(final_X,contrib_H)

#plt.plot(final_X,[0.0]*len(final_X),'go')

#time = np.linspace(0,max_t_val,200)
#uniform_F = [getFt(t,rho,betas) for t in time]
#plt.plot(time,uniform_F)


plt.show()








#plt.plot(x,y)
#plt.plot(x,y,'ro')
#plt.plot(halfway_X,halfway_guess,'go')
#plt.plot(halfway_X,halfway_real,'bo')
#plt.show()





