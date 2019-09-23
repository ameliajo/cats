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
        V0 = 1e-15
        for j in range(1,len(t)):
            sinS = sin(t[j]*AL);
            cosS = cos(t[j]*AL);
            V = (t[j]-t[j-1])*AL;
            if abs(V/V0-1.) > 5e-7:
                if abs(V) <= 0.005:
                  ST = (V*V)/6.-(V**4)/120.0
                  CT = V*0.5-V**3/24.0
                else:
                  sinT = sinS*cosSM-cosS*sinSM;
                  cosT = cosS*cosSM+sinS*sinSM;
                  ST =  1.0-sinT /V;
                  CT = (1.0-cosT)/V;
            EX1 = exp(B*GC[j]);
            EX2 = exp(-A*T*t[j]*t[j]);
            Q2 = cos(B*GS[j])*EX1*EX2-EX2;
            R2 = sin(B*GS[j])*EX1*EX2;
            S[i] +=  Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)-R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM)
            sinSM = sinS;
            cosSM = cosS;
            V0 = V;
            Q1 = Q2;
            R1 = R2;
        S[i] *= F/AL;
    return S;


"""
def SCINT(t,GC,GS,EPS,temp,A,B,F):
    EX1 = exp(B*GC[0])
    EX2 = exp(-A*temp*t[0]**2)
    Q1  = (cos(B*GS[0])*EX1-1.0)*EX2
    R1  = sin(B*GS[0])*EX1*EX2
    S = [0.0]*len(EPS)

    for i in range(len(EPS)):
        AL = A-EPS[i]
        if AL == 0:
            print("AL = 0")
            continue
        S[i] = 0.0
        SM   = t[0]*AL
        SINSM = sin(SM)
        COSSM = cos(SM)
        V0 = 1e-15
        for j in range(1,len(t)):
            U = t[j]*AL
            SINS = sin(U)
            COSS = cos(U)
            V = U-SM
            if abs(V/V0-1.0) > 5e-7:
                if abs(V) <= 0.005:
                    ST = (V**2)/6.0 - (V**2)**2/120.0
                    CT = V*0.5-V**3/24.0
                else:
                    SINT = SINS*COSSM - COSS*SINSM
                    COST = COSS*COSSM + SINS*SINSM
                    ST = 1.0-SINT /V
                    CT =(1.0-COST)/V
            EX1 = exp(B*GC[j])
            EX2 = exp(-A*temp*t[j]**2)
            Q2  = (cos(B*GS[j])*EX1-1.0)*EX2
            R2  = sin(B*GS[j])*EX1*EX2
            S[i] = S[i] + Q2*(ST*SINS+CT*COSS) - Q1*(ST*SINSM-CT*COSSM)
            S[i] = S[i] - R2*(CT*SINS-ST*COSS) + R1*(ST*COSSM+CT*SINSM)
            #if i == 10 and j == 10:
            #    print(S[i])
            SM = U
            SINSM = SINS
            COSSM = COSS
            V0 = V
            Q1 = Q2
            R1 = R2
        S[i] *= F/AL
    #print(S[30:33])
    #print(S[30:33])
    return S




"""


