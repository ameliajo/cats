import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt

def SCINT( t, GC, GS, EPS, T, A, B, F):
    S = [0.0]*len(EPS)

    EX1 = exp(B*GC[0]);
    EX2 = exp(-A*T*t[0]*t[0]);
    Q1 = cos(B*GS[0])*EX1*EX2-EX2;
    R1 = sin(B*GS[0])*EX1*EX2;

    for i in range(len(EPS)):
        AL = A-EPS[i];
        if (AL == 0): 
            print("AL = 0") 
            continue
        S[i] = 0.0
        sinSM = sin(t[0]*AL);
        cosSM = cos(t[0]*AL);
        V0 = 0.0
        invV0 = 1.0e5
        invV  = 0.0
        for j in range(1,len(t)):
            sinS = sin(t[j]*AL);
            cosS = cos(t[j]*AL);
            V = (t[j]-t[j-1])*AL;
            invV = 1.0/V
            if abs(V*invV0-1.) > 5e-7:
                if abs(V) <= 0.005:
                  ST = (V*V)/6.-(V**4)/120.0
                  CT = V*0.5-V**3/24.0
                else:
                  sinT = sinS*cosSM-cosS*sinSM;
                  cosT = cosS*cosSM+sinS*sinSM;
                  ST =  1.0-sinT *invV
                  CT = (1.0-cosT)*invV
            EX1 = exp(B*GC[j]);
            EX2 = exp(-A*T*t[j]*t[j]);
            Q2 = cos(B*GS[j])*EX1*EX2-EX2;
            R2 = sin(B*GS[j])*EX1*EX2;
            S[i] += Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)-\
                    R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM)
            sinSM = sinS;
            cosSM = cosS;
            V0 = V;
            invV0 = invV
            Q1 = Q2;
            R1 = R2;
        S[i] *= F/AL;
    return S;


