#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>


template <typename Float, typename Range>
auto ftrans( Float T, Range X, Range Q, Range t){
    std::vector<double> PC(t.size(),0.0), PS(t.size(),0.0);
    Float SM, SINSM, COSSM, Z, ZM, S, SINS, COSS, U, SINT, COST, ST, CT, H, HM;
    using std::pow; using std::sinh; using std::cosh;
    
    for (size_t i = 1; i < t.size(); ++i){
        SM = X[0]*t[i];
        SINSM = sin(SM);
        COSSM = cos(SM);
        ZM = X[0]*0.5/T;
        for (size_t j = 1; j < X.size(); ++j){
            S = X[j]*t[i];
            SINS = sin(S);
            COSS = cos(S);
            Z = X[j]*0.5/T;
            if (std::abs(S/SM-1.0) > 5e-7){
                U = S-SM;
                if (U <= 0.005){
                    ST = U*U/6.0 - pow(U,4)/120.0;
                    CT = 0.5*U - pow(U,3)/24.0;
                }
                else {
                    SINT = SINS*COSSM - COSS*SINSM;
                    COST = COSS*COSSM + SINS*SINSM;
                    ST =  1.0-SINT /U;
                    CT = (1.0-COST)/U;
                }
            }
            H  = cosh(Z )/sinh(Z );
            HM = cosh(ZM)/sinh(ZM);
            PC[i] = PC[i] + Q[j]/X[j]*H*(ST*SINS+CT*COSS) - Q[j-1]/X[j-1]*HM*(ST*SINSM-CT*COSSM);
            PS[i] = PS[i] + Q[j]/X[j]*  (CT*SINS-ST*COSS) + Q[j-1]/X[j-1]*   (CT*SINSM+ST*COSSM);
            SM = S;
            SINSM = SINS;
            COSSM = COSS;
            ZM = Z;
        }
        PC[i] = PC[i]/t[i];
        PS[i] = PS[i]/t[i];
    }
    return std::make_tuple(PC,PS);




}

