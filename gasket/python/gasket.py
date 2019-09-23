import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
from gtg   import *
from scint import *
from bessl import *
from rconv import *
from acon2 import *

def gasket(X, Q, t, alphas, betas, T, continWgt, freeGasWgt, oscWgt, \
           oscEnergies, oscWgts):
    AM = 1.0/1.0086654; # Convert mass to neutron mass unit
    tbar, GC, GS = GTG( continWgt, T, AM, X, Q, t )
    S1 = [0.0]*len(betas)
    S2 = [0.0]*len(betas)
    S  = [0.0]*len(betas)
    ANK = [0.0]*(len(oscEnergies)*10)
    ARG1 = [0.0,0.0]
    ARG2 = [0.0,0.0]
    BF = [0.0]*10
    NMAX = [0.0]*2
    EPS = [beta*T for beta in betas]
    sab = [0.0]*(len(alphas)*len(betas))

    for a in range(len(alphas)):
        PSQ = alphas[a]*AM*T
        DBW = exp(-PSQ*GC[0])
        APS = PSQ*freeGasWgt/AM
        BPS = PSQ
        DBWP = DBW/3.1416

        S1 = SCINT(t,GC,GS,EPS,T,APS,BPS,DBWP);

        SZCON = DBW*sqrt(AM/(12.566371*PSQ*freeGasWgt*T));
        for i in range(len(betas)):
            S1[i]= S1[i]/exp(betas[i]/2);
            S2[i] = SZCON*exp(-AM*((EPS[i]**2) + (PSQ*freeGasWgt/AM)**2)/(4.*PSQ*T*freeGasWgt));
            S[i] = S1[i]+S2[i];

        for i in range(len(oscEnergies)):
            RR = 0.5*oscEnergies[i]/T;
            U  = exp(RR);
            U  = 0.5*(U-1.0/U);
            ARG1[i] = oscWgt*oscWgts[i]/(AM*oscEnergies[i]*U);
            ARG2[i] = oscWgt*oscWgts[i]/(AM*oscEnergies[i]*sinh(RR)/cosh(RR));
            NMAX[i] = bessl(ARG1[i]*PSQ,BF,10);
            EX = exp(-PSQ*ARG2[i]);
            for j in range(10):
                ANK[i+j*len(oscEnergies)] = BF[j]*EX;
        rconv( NMAX, oscEnergies, ANK, T, S1, betas );
        acon2( NMAX, oscEnergies, ANK, T, SZCON, EPS, AM, freeGasWgt, PSQ, S2 );
        for i in range(len(betas)):
            S[i] = S1[i] + S2[i];
            sab[i+a*len(betas)] = S[i];
    return sab;
