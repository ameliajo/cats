import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt

def ftrans( T, X, Q, t):

    PC = [0.0]*len(t)
    PS = [0.0]*len(t)
    for i in range(1,len(t)):
        sinSM = sin(X[0]*t[i])
        cosSM = cos(X[0]*t[i])
        ZM = X[0]*0.5/T
        for j in range(1,len(X)):
            sinS = sin(X[j]*t[i])
            cosS = cos(X[j]*t[i])
            Z = X[j]*0.5/T
            if abs(X[j]/X[j-1] - 1.0) > 5e-7:
                U = (X[j]*t[i])-(X[j-1]*t[i]);
                if U <= 0.005:
                    ST = U*U/6.0 - U**4/120.0;
                    CT = 0.5*U - U**3/24.0;
                else:
                    SINT = sinS*cosSM - cosS*sinSM;
                    COST = cosS*cosSM + sinS*sinSM;
                    ST =  1.0-SINT /U;
                    CT = (1.0-COST)/U;
            PC[i] += Q[j  ]/X[j  ] * (ST*sinS +CT*cosS ) * cosh(Z) /sinh(Z) - \
                     Q[j-1]/X[j-1] * (ST*sinSM-CT*cosSM) * cosh(ZM)/sinh(ZM);
            PS[i] += Q[j  ]/X[j  ] * (CT*sinS -ST*cosS ) + \
                     Q[j-1]/X[j-1] * (CT*sinSM+ST*cosSM);
            sinSM = sinS;
            cosSM = cosS;
            ZM = Z;
        PC[i] /= t[i];
        PS[i] /= t[i];
    return PC,PS

