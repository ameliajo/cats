#include <iostream>
#include <boost/math/quadrature/gauss.hpp>
#include "getLambda_s.h"

#ifndef CATS_CONTIN_GET_F
#define CATS_CONTIN_GET_F
template <typename Range, typename Float>
auto getF(Range rho, Range betas, Float t, bool useTrapz=false){
  using std::sinh; using std::cosh; using std::exp;
  using namespace boost::math::quadrature;

  // Do the small piece (0-->beta1)
  Float c = rho[1]/(betas[1]*betas[1]);
  Float smallPart = 2.0*c*betas[1] + (c/18.0 - c*t*t/3.0)*pow(betas[1],3)-c*t*t*pow(betas[1],5)/60.0;


  // Do the large piece (beta1-->betaMax)

  std::vector<double> lastBetas (betas.size()-1);
  for (size_t i = 0; i < lastBetas.size(); ++i){lastBetas[i] = betas[i+1];}

  auto func1 = [betas,rho,t](const double& beta) { 
    auto rhoVal = interpolate(betas,rho,beta);
    return rhoVal*cos(beta*t)*(cosh(beta*0.5)/sinh(beta*0.5))/beta; };

  auto func2 = [lastBetas,rho,t](const int& i) { 
    auto beta = lastBetas[i];
    return rho[i]*cos(beta*t)*(cosh(beta*0.5)/sinh(beta*0.5))/beta; };

  Float largePart = useTrapz ? 
      trapz(lastBetas,func2)
    : gauss<double,10>::integrate(func1,betas[1],betas[betas.size()-1]);
  return smallPart + largePart;

}

template <typename Range, typename Float>
auto getG(Range rho, Range betas, Float t, bool useTrapz=false){
  using std::sinh; using std::cosh; using std::exp;
  using namespace boost::math::quadrature;

  // Do the small piece (0-->beta1)
  Float c = rho[1]/(betas[1]*betas[1]);
  Float smallPart = -c/3.0*pow(betas[1],3);


  // Do the large piece (beta1-->betaMax)

  auto func1 = [betas,rho,t](const double& beta) { 
    auto rhoVal = interpolate(betas,rho,beta);
    return rhoVal/beta*sin(beta*t);
  };

  std::vector<double> lastBetas (betas.size()-1);
  for (size_t i = 0; i < lastBetas.size(); ++i){
    lastBetas[i] = betas[i+1];
  }
  auto func2 = [lastBetas,rho,t](const int& i) { 
    auto beta = lastBetas[i];
    return rho[i]/beta*sin(beta*t);
  };

  Float largePart = useTrapz ? -trapz(lastBetas,func2)
                             : -gauss<double,10>::integrate(func1,betas[1],betas[betas.size()-1]);
  return smallPart + largePart;

}

#endif
