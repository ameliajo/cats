import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt,log
from sterp import *


def rconv( nMax, X5, ANK, T, S1, betas):

    SLOG = [0.0]*1000
    SK   = [0.0]*1000
  
    for k in range(len(X5)):
        for i in range(len(S1)):
            if S1[i] <= 0:
                SLOG[i]=-100
            else:
                SLOG[i] = log(S1[i])
        for i in range(len(S1)):
            SK[i]=0.0;
            for j in range(nMax[k]*2+1):
                ind = j-nMax[k];
                BETAIN = abs(betas[i]-float(ind)*X5[k]/T);
                SINT = sterp(BETAIN,betas,SLOG);
                SK[i] += ANK[k+(abs(ind))*len(X5)] * SINT;
            S1[i]=SK[i];




