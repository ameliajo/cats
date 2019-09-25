import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
import matplotlib.pyplot as plt

def ftrans( T, X, Q, time, H):

    PC = [0.0]*len(time)
    PS = [0.0]*len(time)
    for t in range(1,len(time)):
        sinSM = sin(X[0]*time[t])
        cosSM = cos(X[0]*time[t])
        for j in range(1,len(X)):
            sinS = sin(X[j]*time[t])
            cosS = cos(X[j]*time[t])
            if abs(X[j]/X[j-1] - 1.0) > 5e-7:
                U = (X[j]*time[t])-(X[j-1]*time[t]);
                if U <= 0.005:
                    ST = U*U/6.0 - U**4/120.0;
                    CT = 0.5*U - U**3/24.0;
                else:
                    SINT = sinS*cosSM - cosS*sinSM;
                    COST = cosS*cosSM + sinS*sinSM;
                    ST =  1.0-SINT /U;
                    CT = (1.0-COST)/U;
            PC[t] += Q[j  ]/X[j  ] * (ST*sinS +CT*cosS ) * H[j] - \
                     Q[j-1]/X[j-1] * (ST*sinSM-CT*cosSM) * H[j-1];

            PS[t] += Q[j  ]/X[j  ] * (CT*sinS -ST*cosS ) + \
                     Q[j-1]/X[j-1] * (CT*sinSM+ST*cosSM);
            sinSM = sinS;
            cosSM = cosS;
        PC[t] /= time[t];
        PS[t] /= time[t];
    return PC,PS

