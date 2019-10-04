import numpy as np
from numpy import sin,cos,exp,sinh,cosh,sqrt
from gtg   import *
from scint import *
from bessl import *
from rconv import *
from acon2 import *

def gasket(rhoX, rhoY, time, alphas, betas, T, continWgt, freeGasWgt, oscWgt, oscEnergies, oscWgts):
    oscWgts = [oscWgt*oscWgtsVal for oscWgtsVal in oscWgts]
    AM = 1.0/1.0086654; # Convert mass to neutron mass unit
    invT = 1.0/T
    rhoBetas = [rhoX_val*invT for rhoX_val in rhoX]
    print("num timesteps = ",len(time),".    dt = ",time[1]-time[0])
    tbar, H, F = GTG( continWgt, T, AM, rhoX, rhoY, time )
    S2, S = [0.0]*len(betas), [0.0]*len(betas)
    ANK = [0.0]*(len(oscEnergies)*10)
    ARG1, ARG2 = [0.0,0.0], [0.0,0.0]
    ARG2 = [0.0,0.0]
    BF = [0.0]*10
    NMAX = [0.0]*2
    epsilon = [beta*T for beta in betas]
    sab = [0.0]*(len(alphas)*len(betas))

    for a in range(len(alphas)):
        print("doing alpha # ",a," out of ",len(alphas))
        PSQ = alphas[a]*AM*T
        DBW = exp(-PSQ*H[0])

        APS = alphas[a]*T*freeGasWgt
        BPS = alphas[a]*T*AM

        S1 = SCINT(time,H,F,epsilon,T,APS,BPS,DBW);
        print("Got",len(list(filter(lambda x: x < 0, S1)))," negative values")

        SZCON = DBW*sqrt(AM/(12.566371*PSQ*freeGasWgt*T));
        for i in range(len(betas)):
            S1[i]*= exp(-betas[i]*0.5);
            S2[i] = SZCON*exp(-AM*((epsilon[i]**2) + \
                    (PSQ*freeGasWgt/AM)**2)/(4.*PSQ*T*freeGasWgt));

        for i in range(len(oscEnergies)):
            RR = 0.5*oscEnergies[i]*invT;
            U  = exp(RR);
            U  = 0.5*(U-1.0/U);
            ARG1[i] = oscWgts[i]/(AM*oscEnergies[i]*U);
            ARG2[i] = oscWgts[i]/(AM*oscEnergies[i]*sinh(RR)/cosh(RR));
            NMAX[i] = bessl(ARG1[i]*PSQ,BF,10);
            EX = exp(-PSQ*ARG2[i]);
            for j in range(10):
                ANK[i+j*len(oscEnergies)] = BF[j]*EX;
        rconv( NMAX, oscEnergies, ANK, T, S1, betas );
        acon2( NMAX, oscEnergies, ANK, T, SZCON, epsilon, AM, freeGasWgt, PSQ, S2 );
        for i in range(len(betas)):
            #sab[i+a*len(betas)] = S1[i]
            sab[i+a*len(betas)] = S1[i] + S2[i]
    return sab, H, F
