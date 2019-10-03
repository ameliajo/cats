import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
import matplotlib.pyplot as plt

def ftrans( T, X, Q, time, coth):
    H, F = [0.0]*len(time), [0.0]*len(time)

    for t in range(1,len(time)):
        sinSM = sin(X[0]*time[t])
        cosSM = cos(X[0]*time[t])
        for j in range(1,len(X)):
            sinS = sin(X[j]*time[t])
            cosS = cos(X[j]*time[t])
            if abs(X[j]/X[j-1] - 1.0) > 5e-7:
                theta = (X[j]*time[t])-(X[j-1]*time[t]);
                if theta <= 0.005:
                    ST = theta*theta/6.0 - theta**4/120.0;
                    CT = 0.5*theta - theta**3/24.0;
                else:
                    SINT = sinS*cosSM - cosS*sinSM;
                    COST = cosS*cosSM + sinS*sinSM;
                    ST =  1.0-SINT /theta;
                    CT = (1.0-COST)/theta;
            H[t] += Q[j  ]/X[j  ] * (ST*sinS +CT*cosS ) * coth[j] - \
                    Q[j-1]/X[j-1] * (ST*sinSM-CT*cosSM) * coth[j-1];

            F[t] += Q[j  ]/X[j  ] * (CT*sinS -ST*cosS ) + \
                    Q[j-1]/X[j-1] * (CT*sinSM+ST*cosSM);
            sinSM = sinS;
            cosSM = cosS;
        H[t] /= time[t];
        F[t] /= time[t];
    return H,F

