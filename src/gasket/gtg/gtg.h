#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>
#include "helpers_gtg/intg.h"
#include "helpers_gtg/ftrans.h"

template <typename Range, typename Float>
auto GTG(Float wgt, Float T, Float AM, Range X, Range Q, Range t ){
    using std::pow; using std::cosh; using std::sinh;

//SUBROUTINE GTG(wgt,T,AM,X,Q,t,PC,PS,JS3,NT,TBAR)

//INTEGER :: JS3, NT
//REAL(8) :: X(JS3), Q(JS3), t(NT)
//REAL(8) :: X2(5), Q2(5), t2(10)
//REAL(8) :: wgt, T, AM, TBAR

    Float F, H, C, CS, S, invU, TBAR;

    Float U = Q[0]*X[0]/3.0;
    // puts integral result in variable A -- JS3 is the length of X and Q
    //INTG(X,Q,A,JS3);
    Float A = INTG(X,Q);

    Float norm = wgt / AM / (U+A);
    auto out = ftrans(T, X, Q, t);
    Range PC = std::get<0>(out);
    Range PS = std::get<1>(out);

    F = X[0]*0.5/T;
    H = cosh(F)/sinh(F);
    for (size_t i = 1; i < t.size(); ++i ){
      U = X[0]*t[i];
      if (U <= 0.005) {
        C  = 0.5*U - pow(U,3)/24.0;
        S  = U - pow(U,3)/6.0;
        CS = U/3.0 - pow(U,3)/30.0;
      } else {
        C = cos(U);
        S = sin(U);
        invU = 1.0/U;
        CS = S*invU*invU - C*invU;
        C = (1.-C)*invU;
      }
      PC[i]=PC[i]+Q[0]/U*(H*(S-C)+C/F);
      PS[i]=PS[i]+Q[0]*CS;
    }
    // Putting in the t=0 terms
    PC[0] = 0.5*Q[0]*(1./F+H);
    PS[0] = 0.0;
    TBAR = Q[0]*T*X[0]/3.0;

    for (size_t i=0; i < X.size(); ++i ){
      F = X[i]*0.5/T;
      H = cosh(F)/sinh(F);
      Q[i] = Q[i]*H/X[i];
    }
    A = INTG(X,Q);
    PC[0]=PC[0]+A;
    for (size_t i = 0; i < X.size(); ++i ){
      Q[i] = Q[i]*pow(X[i],2)*0.5;
    }
    A = INTG(X,Q);
    TBAR += A;

    for (size_t i = 0; i < t.size(); ++i ){
      PC[i] = PC[i]*norm;
      PS[i] = PS[i]*norm;
    }
    TBAR = TBAR*norm*AM/wgt;
    return std::make_tuple(TBAR,PC,PS);
}


