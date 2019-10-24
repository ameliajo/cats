import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt

def sterp( B, betas, sLog ):
    if B <= betas[0]:
        return exp(sLog[0])
    elif B >= betas[len(betas)-1]:
        return 0.0
    for IC in range(len(betas)-1):
        if B >= betas[IC] and B <= betas[IC+1]:
            SL = sLog[IC]+((B-betas[IC])/(betas[IC+1]-betas[IC]))*(sLog[IC+1]-sLog[IC]);
            return exp(SL);
    return 0.0





