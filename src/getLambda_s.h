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

  normalizeRho(betas, rho, useOld);

  auto f1 = [betas,rho](int i) { 
    auto P = (abs(betas[i]) < 1e-12) ? rho[1]/(betas[1]*betas[1]) 
                                     : rho[i]/(2.0*betas[i]*sinh(betas[i]*0.5));
    return P*2.0*cosh(betas[i]*0.5);
  };
  auto f2 = [betas,rho](const double& beta) { 
    auto rhoVal = interpolate(betas,rho,std::abs(beta));
    auto P = (abs(beta) < 1e-12) ? rho[1]/(betas[1]*betas[1]) 
                                 : rhoVal/(2.0*beta*sinh(beta*0.5));
    return P*2.0*cosh(beta*0.5);
  };

  std::cout << trapz(betas,f1)<<std::endl;
  std::cout << gauss<double,100>::integrate(f2,0.0,betas[betas.size()-1]) <<std::endl;
  return useOld ? trapz(betas,f1)
                : gauss<double,100>::integrate(f2,0.0,betas[betas.size()-1]);

}

#endif 



