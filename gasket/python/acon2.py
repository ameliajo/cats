import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt

def acon2(NMAX, X5, ANK, T, SZCON, EPS, AM, W1, PSQ, S2):
    SK = [0.0]*1000
    SINT_denominator = AM/(4.0*PSQ*T*W1)
    PSQ_power = (PSQ*W1/AM)**2
    for i in range(len(S2)):
        for j in range(NMAX[1]*2+1):
            for k in range(NMAX[0]*2+1):
                EIN = abs(EPS[i] - (k-NMAX[0])*X5[0]\
                                 - (j-NMAX[1])*X5[1])
                SINT = SZCON * exp(-(EIN**2+PSQ_power)*SINT_denominator)
                SK[i] += ANK[1+(abs(j-NMAX[1]))*len(X5)] * \
                         ANK[0+(abs(k-NMAX[0]))*len(X5)] * SINT;
        S2[i] = SK[i]

