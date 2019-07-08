#include <iostream>
#include <boost/math/quadrature/gauss.hpp>

template <typename Range, typename Float>
auto interpolate(Range betas, Range rho, Float t){
  for (size_t i = 0; i < betas.size()-1; ++i){
    if (betas[i] <= t and t <= betas[i+1]){
      return (rho[i+1]-rho[i])/(betas[i+1]-betas[i])*(t-betas[i])+rho[i];
    }
  }
  return 0.0;
}

template <typename Range,typename Function>
auto trapz(Range x, Function func){
  double integral = 0.0;
  for (size_t i = 0; i < x.size()-1; ++i){
    integral += 0.5*(func(i)+func(i+1))*(x[i+1]-x[i]);
  }
  return integral;
}

template <typename Range>
void normalizeRho(const Range& betas, Range& rho, bool useOld){
  using namespace boost::math::quadrature;
  auto f1 = [betas,rho](auto i) { return rho[i]; };
  auto f2 = [betas,rho](const double& beta) { return interpolate(betas,rho,beta); };
  double invArea = useOld ? 
      1.0/(trapz(betas,f1))//-0.5*rho[1]*betas[1] + rho[1]*pow(betas[1],3)/3.0)
   :  1.0/(gauss<double,100>::integrate(f2,0.0,betas[betas.size()-1]));
  for (auto& x : rho){ x *= invArea; }
}




template <typename Range>
auto getLambda_s(Range betas, Range rho, bool useOld = false){
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

  return useOld ? trapz(betas,f1)
                : gauss<double,10>::integrate(f2,0.0,betas[betas.size()-1]);

}





