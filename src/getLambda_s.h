#ifndef CATS_GET_LAMBDA_S
#define CATS_GET_LAMBDA_S
#include <iostream>
#include <boost/math/quadrature/gauss.hpp>
#include "generalTools/trapezoid.h"
#include "generalTools/interpolate.h"
#include "generalTools/normalize.h"

template <typename Range>
auto getLambda_s(Range betas, Range rho, bool useOld=false){
  using std::sinh; using std::cosh; using std::exp;
  using namespace boost::math::quadrature;

  auto maxBeta = betas[betas.size()-1];
  normalizeRho(betas, rho, useOld);

  std::vector<double> P(rho.size(),0.0);
  P[0] = rho[1]/(betas[1]*betas[1]);
  for ( size_t i = 1; i < P.size(); ++i ){ 
    P[i] = rho[i]/(2.0*betas[i]*sinh(betas[i]*0.5));
  }
  //std::cout << P[0] << "  " << P[1] << "  " << P[2] << std::endl;
  //std::cout << P[3] << "  " << P[4] << "  " << P[5] << std::endl;


  auto trapzIntegrand   = [betas,P](int i){return P[i]*2.0*cosh(betas[i]*0.5);};
  auto gaussLegendreArg = [betas,P](auto beta){ 
    return exp(-beta*0.5)*interpolate(betas,P,abs(beta)); 
  };

  //std::cout << gauss<double,10>::integrate(gaussLegendreArg,-maxBeta,maxBeta) << std::endl;
  return useOld ? trapz(betas,trapzIntegrand)
                : gauss<double,10>::integrate(gaussLegendreArg,-maxBeta,maxBeta);
}

#endif 



