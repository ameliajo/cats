import numpy as np
from numpy import sinh, cosh, exp



def getConvolAtPoint(i, delta, t1, t2):
    sumVal = 0.0
    for j in range(-len(t1)+1,len(t1)):
        if ( i - j >= len(t2) ):
            return sumVal
        expVal = exp(   j *delta) if j < 0 else \
                 exp((i-j)*delta) if i < j else 1.0
        toAdd = t1[abs(j)]*t2[abs(i-j)]*expVal
        sumVal += 0.5*toAdd if ( j == -len(t1)+1 or j == len(t1)-1 ) else toAdd
    return sumVal

def convol(t1, t2, delta, nn):
    t3 = [0.0]*nn
    for i in range(nn):
        t3[i] = getConvolAtPoint(i,delta,t1,t2)*delta
        if t3[i] < 1e-30:
            t3[i] = 0.0
    return t3


