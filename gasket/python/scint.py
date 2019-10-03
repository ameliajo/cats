import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
import time
from math import pi

def SCINT( t, H, F, epsilon, T, w1_alpha, alpha, DBW):
    S = [0.0]*len(epsilon)

    EX1 = exp(alpha*H[0]);
    EX2 = exp(-w1_alpha*T*t[0]*t[0]);
    Q1 = cos(alpha*F[0])*EX1*EX2-EX2;
    R1 = sin(alpha*F[0])*EX1*EX2;

    #start= time.time()
    for i in range(len(epsilon)):
        AL = w1_alpha-epsilon[i];
        #print(A,AL)
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
            EX1 = exp(alpha*H[j]);
            EX2 = exp(-w1_alpha*T*t[j]*t[j]);
            Q2 = cos(alpha*F[j])*EX1*EX2-EX2;
            #print('{:.2e}'.format(Q1)>'{:.2e}'.format(Q2))
            R2 = sin(alpha*F[j])*EX1*EX2;
            #print( ( Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)-\
            #         R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM) ) /\
            #     abs(Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)-\
            #         R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM) ) )

            S[i] += Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)-\
                    R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM)
            sinSM = sinS;
            cosSM = cosS;
            V0 = V;
            Q1 = Q2;
            R1 = R2;
        S[i] *= DBW/(pi*AL);
        #print(A,AL,S[i]/abs(S[i]))
    #end = time.time()
    #print("-----",end-start)
    return S;


