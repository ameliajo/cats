#include "continTools/get_F_G.h"
#include "getLambda_s.h"
#include <boost/math/quadrature/gauss.hpp>
#include "generalTools/gaussLaguerre.h"

template <typename Range, typename Float>
auto getSab(Range betas, Range rho, Float alpha, Float beta){
  using std::exp; using std::cos;
  auto lambda_s = getLambda_s(betas,rho,false);
  auto constFactors = exp(-alpha*lambda_s)/M_PI;
  std::cout << std::endl;
  std::cout << constFactors << std::endl;
  std::cout << "lambda_s = " << lambda_s << std::endl;
  std::cout << std::endl;

  //auto integrand = [rho,betas,alpha,beta](auto t){
  //  Float F = getF(rho,betas,t,true);
  //  Float G = getG(rho,betas,t,true);
  //  return exp(alpha*F)*cos(beta*t-alpha*G);
  //};
  //std::cout << boost::math::quadrature::gauss<double,10>::integrate(integrand,0.0,10.0) << std::endl;
  
  int N = 10;
  auto points  = std::get<0>(gaussLaguerre(N,0.0,2.0));
  auto weights = std::get<1>(gaussLaguerre(N,0.0,2.0));
  double contribs = 0.0;
  for (int i = 0; i < N; ++i){
    auto t = points[i];
    auto F = getF(rho,betas,t,false);
    auto G = getG(rho,betas,t,false);
    contribs += weights[i]*constFactors*exp(alpha*F+t)*cos(beta*t-alpha*G);
  }
  return contribs;

}



