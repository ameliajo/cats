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
  //for ( auto x : rho){ std::cout << x << " ";}
  //std::cout << std::endl;

  auto f0 = [betas,rho](int i) { 
    return i == 0 ? rho[1]/(betas[1]*betas[1]) 
                  : rho[i]/(2.0*betas[i]*sinh(betas[i]*0.5));
  };

  auto f1 = [betas,rho](int i) { 
    auto P = (abs(betas[i]) < 1e-12) ? rho[1]/(betas[1]*betas[1]) 
                                     : rho[i]/(2.0*betas[i]*sinh(betas[i]*0.5));
    return P*2.0*cosh(betas[i]*0.5);
  };
  auto f2 = [betas,rho](const double& beta) { 
    auto rhoVal = interpolate(betas,rho,std::abs(beta));
    
    auto P = (abs(beta) < 1e-12) ? rho[1]/(betas[1]*betas[1]) 
                                 : rhoVal/(2.0*beta*sinh(beta*0.5));
    return P*exp(-beta*0.5);
  };

  std::vector<double> P_vec(rho.size(),0.0);
  for ( size_t i = 0; i < P_vec.size(); ++i ){ P_vec[i] = f0(i); }
  auto f3 = [betas,P_vec](auto beta){ 
    //std::cout << "-------        "<<interpolate(betas,P_vec,abs(beta))<<std::endl;
    return exp(-beta*0.5)*betas[betas.size()-1]*interpolate(betas,P_vec,abs(beta)); };

  //std::cout << trapz(betas,f1)<<std::endl;
  //std::cout << gauss<double,10>::integrate(f2,-betas[betas.size()-1],betas[betas.size()-1]) <<std::endl;
  //std::cout << "-----   "<<P_vec[0] << std::endl;
  //std::cout << "-----   "<<P_vec[1] << std::endl;
  //std::cout << "-----   "<<P_vec[2] << std::endl;
  //std::cout << "-----   "<<P_vec[3] << std::endl;
  //std::cout << "-----   "<<P_vec[4] << std::endl;
  //std::cout << std::endl;
  return useOld ? trapz(betas,f1)
                : gauss<double,10>::integrate(f3,-betas[betas.size()-1],betas[betas.size()-1]);
                //: gauss<double,10>::integrate(f2,-betas[betas.size()-1],betas[betas.size()-1]);

}

#endif 



