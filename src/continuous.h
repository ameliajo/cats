#include "continTools/get_F_G.h"
#include "getLambda_s.h"
#include <boost/math/quadrature/gauss.hpp>



template <typename Range, typename Float>
auto getSab(Range betas, Range rho, Float alpha, Float beta){
  using std::exp; using std::cos;
  auto lambda_s = getLambda_s(betas,rho,true);
  auto sab = exp(-alpha*lambda_s)/M_PI;

  auto integrand = [rho,betas,alpha,beta](auto t){
    Float F = getF(rho,betas,t,true);
    Float G = getG(rho,betas,t,true);
    return exp(alpha*F)*cos(beta*t-alpha*G);
  };

  auto f1 = [](auto x){return 3*x;};
  std::cout << boost::math::quadrature::gauss<double,10>::integrate(f1,0.0,10.0) << std::endl;

  return;
  
  std::cout << rho[0]<<betas[0]<<alpha<<beta<<std::endl;

}



