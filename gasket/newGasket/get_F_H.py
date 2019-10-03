import numpy as np
from numpy import sin,cos,sinh,cosh
import matplotlib.pyplot as plt
from waterDataContinuous import *

def get_CT(x):
    return (1.0-cos(x))/x if x > 0.005 else x*0.5 - x**3/24.0
def get_ST(x):
    return  1.0-sin(x) /x if x > 0.005 else x**2/6.0 - x**4/120.0
def get_sin(x):
    return sin(x) if x > 0.005 else x - x**3*0.1666666
def get_F0(x):
    return (sin(x) - x*cos(x))/x**2 if x > 0.005 else x*0.3333333 - x**3*0.0333333


def get_F_H(betas,Q,time,experiment=False):
    coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in betas]
    F, H = [0.0]*len(time), [0.0]*len(time)

    U = Q[0]*betas[0]/3.0
    A = np.trapz(Q,betas)

    norm = 1.0 / (U+A)

    for t in range(1,len(time)):
        #print(t/len(time)*100)
        sinNext = sin(betas[0]*time[t])
        cosNext = cos(betas[0]*time[t])

        for j in range(1,len(betas)):
            sinLast = sinNext; sinNext = sin(betas[j]*time[t]) # Advance sin and cos 
            cosLast = cosNext; cosNext = cos(betas[j]*time[t])

            if abs(betas[j]/betas[j-1] - 1.0) > 5e-7:
                theta = (betas[j]-betas[j-1])*time[t];
                ST, CT = get_ST(theta), get_CT(theta)

            H[t] += Q[j  ]/betas[j  ] * (ST*sinNext + CT*cosNext ) * coth[j  ] - \
                    Q[j-1]/betas[j-1] * (ST*sinLast - CT*cosLast ) * coth[j-1]
            F[t] += Q[j  ]/betas[j  ] * (CT*sinNext - ST*cosNext ) + \
                    Q[j-1]/betas[j-1] * (CT*sinLast + ST*cosLast )

        H[t] /= time[t];
        F[t] /= time[t];
 
    """
    plt.plot(time,F,label='F')
    plt.plot(time,H,label='H')
    plt.legend(loc='best')
    plt.show()
    H2 = H[:]
    F2 = F[:]
    """

    for i in range(1,len(time)):
        U = betas[0]*time[i]
        if U <= 0.005:
            #C  = 0.5*U - (U**3)/24.0;
            #S  = U - (U**3)/6.0;
            CS = U/3.0 - (U**3)/30.0;
        else:
            #C = cos(U);
            #S = sin(U);
            invU = 1.0/U;
            CS = sin(U)*invU*invU - cos(U)*invU;
            #CS = S*invU*invU - C*invU;
            #C = (1.-C)*invU;
        C = get_CT(betas[0]*time[i])
        S = get_sin(betas[0]*time[i])
        #CS = get_F0(betas[0]*time[i])
        #if abs(S2-S)>1e-6:
        #    print("S :( ")
        #    exit()
        #if abs(C2-C)>1e-6:
        #    print("C :( ")
        #    exit()
        #if abs(CS2-CS)>1e-6:
        #    print(U)
        #    print("CS :( ")
        #    exit()
        H[i] += Q[0]/U*(coth[0]*(S-C)+C/(betas[0]*0.5));
        F[i] += Q[0]*CS;


    # Putting in the t=0 terms
    H[0] = Q[0]/betas[0] + 0.5*Q[0]*coth[0]
    F[0] = 0.0;

    """
    for i in range(len(betas)):
        Q[i] = Q[i]*coth[i]/betas[i];
    A = np.trapz(Q,betas)
    H[0] += A
    for i in range(len(betas)):
      Q[i] = Q[i]*(betas[i]**2)*0.5;

    A = np.trapz(Q,betas)

    for i in range(len(time)):
      H[i] = H[i]*norm;
      F[i] = F[i]*norm;
    """

    H[0] += np.trapz([Q[i]*coth[i]/betas[i] for i in range(len(betas))],betas)
    for i in range(len(betas)):
      Q[i] = Q[i]*(betas[i]**2)*0.5;

    for i in range(len(time)):
      H[i] *= norm;
      F[i] *= norm;


    return F,H


if __name__=="__main__":
    invT = 1.0/0.025
    time = np.linspace(0,100,1e4)
    betas = [rhoX_val*invT for rhoX_val in X]
    
    F,H = get_F_H(betas,Q,time)
    
    #invArea = 1.0/np.trapz(F,x=time); F = [x*invArea for x in F]
    #invArea = 1.0/np.trapz(H,x=time); H = [x*invArea for x in H]
 
    
    plt.plot(time,H,label='H')
    plt.plot(time,F,label='F')
    plt.legend(loc='best')
    plt.show()



