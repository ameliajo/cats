import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt,log

def bessl( xVal, B, NX ):
    BF = [0.0]*60
    NMAX = 0;

    if xVal < 0.05:
        FNFACT = 1.0
        X2N    = 1.0
        X2     =xVal*0.5
        for i in range(NX):
            B[i] = X2N/FNFACT;
            X2N *= X2;
            FNFACT = FNFACT*float(i+1);
            if B[i] < 1e-20:
                B[i] = 0.0;
                if (NMAX == 0):
                    NMAX = i
    else:
        for i in range(NX):
            B[i]=0.0;
            if xVal-1.0 < 0:
                IORD=int(-37.0/(.43429*log(.1*xVal)))
                if IORD-5 < 0:
                    B[0]=1
            else:
                IORD=30;
            BF[IORD-2] = 1.e-37;
            BF[IORD-1] = 0.0;
            TA=1.e-37;
            for j in range(2,IORD):
                BF[IORD-j-1]=float(IORD-j)*2.0*BF[IORD-j]/xVal+BF[IORD-j+1];
                TA=TA+BF[IORD-j-1];
            TA=2.0*TA-BF[0];
            FACT=exp(xVal)/TA;
            for j in range(NX):
                B[j]=FACT*BF[j];
                if B[j] < 1e-20:
                    B[j] = 0.0;
                    if NMAX == 0:
                        NMAX = j-1+1; 
    return NMAX;
