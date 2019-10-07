import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,exp


def func(x):
    return 5.0*x*exp(-x)+exp(-x)
def get_CT(x):
    return (1.0-cos(x))/x if x > 0.005 else x*0.5 - x**3/24.0
def get_ST(x):
    return  1.0-sin(x) /x if x > 0.005 else x**2/6.0 - x**4/120.0
def get_sin(x):
    return sin(x) if x > 0.005 else x - x**3*0.1666666
def get_F0(x):
    return sin(x)/x**2 - cos(x)/x if x > 0.005 else x/3.0 - (x**3)/30.0
 

rho = [func(x) for x in np.linspace(0,10,50)]
rho = [0.0,rho[0]*0.5] + rho + [0.0]
betas = np.linspace(1e-6,3.0,len(rho))
#plt.plot(betas,rho,'ro'); plt.show(); exit()

def getFs(t,rho,betas):
    if t == 0.0: return 0.0
    Fs = rho[0]*t*get_F0(betas[0]*t)
    sinNext = sin(betas[0]*t)
    cosNext = cos(betas[0]*t)

    for j in range(1,len(betas)):
        sinLast = sinNext; sinNext = sin(betas[j]*t) # Advance sin
        cosLast = cosNext; cosNext = cos(betas[j]*t) # and cos
        if abs(betas[j]/betas[j-1] - 1.0) > 5e-7:
            theta = (betas[j]-betas[j-1])*t;
            ST, CT = get_ST(theta), get_CT(theta)
        Fs += rho[j  ]/betas[j  ] * (CT*sinNext - ST*cosNext ) + \
              rho[j-1]/betas[j-1] * (CT*sinLast + ST*cosLast )
    return Fs/t

def getF(rho,betas,spacing,max_t_val,tolerance):
    xL = 0.0;        yL = getFs(xL,rho,betas); final_X = [xL]
    xR = xL+spacing; yR = getFs(xR,rho,betas); final_Y = [yL]
    while xR < max_t_val:
        xM = (xL+xR)*0.5
        yM = getFs(xM,rho,betas)
        guess = (yL+yR)*0.5
        if abs(guess-yM) < tolerance:
            final_X.append(xR)
            final_Y.append(yR)
            xL = xR; xR = xL + spacing
            yL = yR; yR = getFs(xR,rho,betas)
        else:
            xR = xM; yR = yM
    return final_X,final_Y

spacing = 1.0
max_t_val = 100.0
tolerance = 1e-3
final_X, final_Y = getF(rho,betas,spacing,max_t_val,tolerance)

plt.plot(final_X,final_Y,'ro')
plt.plot(final_X,[0.0]*len(final_X),'go')
plt.show()





#plt.plot(x,y)
#plt.plot(x,y,'ro')
#plt.plot(halfway_X,halfway_guess,'go')
#plt.plot(halfway_X,halfway_real,'bo')
#plt.show()





