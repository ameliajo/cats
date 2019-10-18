#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>


template <typename Range, typename Float>
auto scint(Range t,Range GC,Range GS,Range EPS, Float T, Float A, Float B , Float F){
  Range S (EPS.size(),0.0);
  using std::exp; using std::pow;

  Float SM, sinSM, cosSM, V0, U, sinS, cosS, V, ST, CT, sinT, cosT, AL, EX1, EX2, Q1, Q2, R1, R2;

  EX1 = exp(B*GC[0]);
  EX2 = exp(-A*T*t[0]*t[0]);
  Q1 = cos(B*GS[0])*EX1*EX2-EX2;
  R1 = sin(B*GS[0])*EX1*EX2;

  for (size_t i = 0; i < EPS.size(); ++i){
    AL = A-EPS[i];
    if (AL == 0){ std::cout << "AL = 0" << std::endl; continue; }

    S[i] = 0.0;
    sinSM = sin(t[0]*AL);
    cosSM = cos(t[0]*AL);
    V0 = 0.0;

    for (size_t j = 1; j < t.size(); ++j){
      sinS = sin(t[j]*AL);
      cosS = cos(t[j]*AL);
      V = (t[j]-t[j-1])*AL;
      if (abs(V/V0-1.) > 5e-7){
        if (abs(V) <= 0.005){
          ST = (V*V)/6.-(pow(V,4))/120.;
          CT = V*0.5-pow(V,3)/24.;
        } else {
          sinT = sinS*cosSM-cosS*sinSM;
          cosT = cosS*cosSM+sinS*sinSM;
          ST =  1.0-sinT /V;
          CT = (1.0-cosT)/V;
        }
      }
      EX1 = exp(B*GC[j]);
      EX2 = exp(-A*T*t[j]*t[j]);
      Q2 = cos(B*GS[j])*EX1*EX2-EX2;
      R2 = sin(B*GS[j])*EX1*EX2;
      S[i] +=  Q2*(ST*sinS+CT*cosS)-Q1*(ST*sinSM-CT*cosSM)
              -R2*(CT*sinS-ST*cosS)-R1*(ST*cosSM+CT*sinSM);
      sinSM = sinS;
      cosSM = cosS;
      V0 = V;
      Q1 = Q2;
      R1 = R2;
    }
    S[i] *= F/AL;
  }
  return S;
}

