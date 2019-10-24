import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
from ftrans import *

def GTG(wgt, T, AM, X, Q, t ):
    invT = 1.0/T
    rhoBetas = [rhoX_val*invT for rhoX_val in X]
    coth = [cosh(beta*0.5)/sinh(beta*0.5) for beta in rhoBetas]

    U = Q[0]*X[0]/3.0
    A = np.trapz(Q,X)

    norm = wgt / AM / (U+A)
    H, F = ftrans(T, X, Q, t, coth)

    for i in range(1,len(t)):
        U = X[0]*t[i];
        if U <= 0.005:
            C  = 0.5*U - (U**3)/24.0;
            S  = U - (U**3)/6.0;
            CS = U/3.0 - (U**3)/30.0;
        else:
            C = cos(U);
            S = sin(U);
            invU = 1.0/U;
            CS = S*invU*invU - C*invU;
            C = (1.-C)*invU;
        H[i] += Q[0]/U*(coth[0]*(S-C)+C/(rhoBetas[0]*0.5));
        F[i] += Q[0]*CS;

    # Putting in the t=0 terms
    H[0] = 0.5*Q[0]*(1./(rhoBetas[0]*0.5)+coth[0]);
    F[0] = 0.0;
    TBAR = Q[0]*T*X[0]/3.0;

    for i in range(len(X)):
        Q[i] = Q[i]*coth[i]/X[i];
    A = np.trapz(Q,X)
    H[0]=H[0]+A;
    for i in range(len(X)):
      Q[i] = Q[i]*(X[i]**2)*0.5;

    A = np.trapz(Q,X)
    TBAR += A;

    for i in range(len(t)):
      H[i] = H[i]*norm;
      F[i] = F[i]*norm;

    TBAR = TBAR*norm*AM/wgt;
    return TBAR,H,F


