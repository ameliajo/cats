#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>
#include "helperFiles/gtg.h"
#include "helperFiles/scint.h"
#include "helperFiles/bessl.h"
#include "helperFiles/rconv.h"
#include "helperFiles/acon2.h"

template <typename Range, typename Float>
auto gasket(Range X, Range Q, Range t, Range alphas, Range betas, Float T, 
  Float continWgt, Float freeGasWgt, Float oscWgt, Range oscEnergies, Range oscWgts){
  using std::exp; using std::pow; using std::sinh; using std::cosh;
  Float AM = 1.0/1.0086654; // Convert mass to neutron mass unit
  Float APS, BPS, DBWP, DBW, PSQ, SZCON, RR, U, EX;

  auto out = GTG( continWgt, T, AM, X, Q, t );
  Float tbar = std::get<0>(out);
  Range GC = std::get<1>(out);
  Range GS = std::get<2>(out);

  Range S1(betas.size(),0.0), 
        S2(betas.size(),0.0), 
         S(betas.size(),0.0);
  Range ANK(oscEnergies.size()*10), ARG1(2), ARG2(2), BF(10), NMAX(2);
  Range EPS = betas | ranges::view::transform([T](auto beta){return beta*T;});
  Range sab (alphas.size()*betas.size());


  for (size_t a = 0; a < alphas.size(); ++a){
    PSQ = alphas[a]*AM*T;
    DBW = exp(-PSQ*GC[0]);
    APS = PSQ*freeGasWgt/AM;
    BPS = PSQ;

    DBWP = DBW/3.1416;

    S1 = scint(t,GC,GS,EPS,T,APS,BPS,DBWP);

    SZCON = DBW*sqrt(AM/(12.566371*PSQ*freeGasWgt*T));
    for (size_t i = 0; i < betas.size(); ++i ){
      S1[i]= S1[i]/exp(betas[i]/2);
      S2[i] = SZCON*exp(-AM*(pow(EPS[i],2) + 
              pow(PSQ*freeGasWgt/AM,2))/(4.*PSQ*T*freeGasWgt));
      S[i] = S1[i]+S2[i];
    }

    for (size_t i = 0; i < oscEnergies.size(); ++i ){
      RR = 0.5*oscEnergies[i]/T;
      U  = exp(RR);
      U  = 0.5*(U-1.0/U);
      ARG1[i] = oscWgt*oscWgts[i]/(AM*oscEnergies[i]*U);
      ARG2[i] = oscWgt*oscWgts[i]/(AM*oscEnergies[i]*sinh(RR)/cosh(RR));
      NMAX[i] = bessl(ARG1[i]*PSQ,BF,10);
      EX = exp(-PSQ*ARG2[i]);
      for (size_t j = 0; j < 10; ++j){
        ANK[i+j*oscEnergies.size()] = BF[j]*EX;
      }
    }
    rconv( NMAX, oscEnergies, ANK, T, S1, betas );
    acon2( NMAX, oscEnergies, ANK, T, SZCON, EPS, AM, freeGasWgt, PSQ, S2 );
    for ( size_t i = 0; i < betas.size(); ++i ){
      S[i] = S1[i] + S2[i];
      sab[i+a*betas.size()] = S[i];
    }
  }
  return sab;

}
