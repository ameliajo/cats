import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
from ftrans import *


def GTG(wgt, T, AM, X, Q, t ):

    U = Q[0]*X[0]/3.0
    # puts integral result in variable A -- JS3 is the length of X and Q
    # INTG(X,Q,A,JS3);
    A = np.trapz(Q,X)

    norm = wgt / AM / (U+A)
    PC, PS = ftrans(T, X, Q, t)

    F = X[0]*0.5/T;
    H = cosh(F)/sinh(F);
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
        PC[i]=PC[i]+Q[0]/U*(H*(S-C)+C/F);
        PS[i]=PS[i]+Q[0]*CS;
    # Putting in the t=0 terms
    PC[0] = 0.5*Q[0]*(1./F+H);
    PS[0] = 0.0;
    TBAR = Q[0]*T*X[0]/3.0;

    for i in range(len(X)):
        F = X[i]*0.5/T;
        H = cosh(F)/sinh(F);
        Q[i] = Q[i]*H/X[i];
    A = np.trapz(Q,X)
    PC[0]=PC[0]+A;
    for i in range(len(X)):
      Q[i] = Q[i]*(X[i]**2)*0.5;

    A = np.trapz(Q,X)
    TBAR += A;

    for i in range(len(t)):
      PC[i] = PC[i]*norm;
      PS[i] = PS[i]*norm;

    TBAR = TBAR*norm*AM/wgt;
    return TBAR,PC,PS


