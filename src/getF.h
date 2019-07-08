#include <iostream>
#include <boost/math/quadrature/gauss.hpp>
#include "getLambda_s.h"

template <typename Range, typename Float>
auto getF(Range rho, Range betas, Float beta, Float t){
  using std::sinh; using std::cosh; using std::exp;
  using namespace boost::math::quadrature;

  if ( abs(beta) < 1e-12 ){ 
    auto c = rho[1]/(betas[1]*betas[1]);
    return 2.0*c*betas[1] + (c/18.0 - c*t*t/3.0)*pow(betas[1],3)-c*t*t*pow(betas[1],5)/60.0;
  } 

  auto function = [betas,rho,t](const double& beta) { 
    auto rhoVal = interpolate(betas,rho,beta);
    return rhoVal*cos(beta*t)*(cosh(beta*0.5)/sinh(beta*0.5))/beta;
  };

  return gauss<double,10>::integrate(function,betas[1],betas[betas.size()-1]);





}

