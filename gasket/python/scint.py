import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
import time

def SCINT( t, GC, GS, EPS, T, A, B, DBWP):
    S = [0.0]*len(EPS)

    EX1 = exp(B*GC[0]);
    EX2 = exp(-A*T*t[0]*t[0]);
    Q1 = cos(B*GS[0])*EX1*EX2-EX2;
    R1 = sin(B*GS[0])*EX1*EX2;

    #start= time.time()
    for i in range(len(EPS)):
        AL = A-EPS[i];
        if (AL == 0): print("AL = 0"); continue
        sinSM = sin(t[0]*AL);
        cosSM = cos(t[0]*AL);
        V0 = 1e-5
        for j in range(1,len(t)):
            sinS = sin(t[j]*AL);
            cosS = cos(t[j]*AL);
            V = (t[j]-t[j-1])*AL;
            if abs(V/V0-1.) > 5e-7:
                if abs(V) <= 0.005:
                  ST = (V*V)*0.166666666-(V**4)*0.00833333333
                  CT = V*0.5-V**3*0.0416666666
                else:
                  sinT = sinS*cosSM-cosS*sinSM;
                  cosT = cosS*cosSM+sinS*sinSM;
                  ST =  1.0-sinT /V
                  CT = (1.0-cosT)/V
            EX1 = exp(B*GC[j]);
            EX2 = exp(-A*T*t[j]*t[j]);
            Q2 = cos(B*GS[j])*EX1*EX2-EX2;
            R2 = sin(B*GS[j])*EX1*EX2;
            S[i] += Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)-\
                    R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM)
            sinSM = sinS;
            cosSM = cosS;
            V0 = V;
            Q1 = Q2;
            R1 = R2;
        S[i] *= DBWP/AL;
    #end = time.time()
    #print("-----",end-start)
    return S;


