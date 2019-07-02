#include <iostream>
#include <boost/math/quadrature/gauss.hpp>


template <typename Range, typename Float>
auto getValue(Range betas, Range rho, Float t){
  for (size_t i = 0; i < betas.size()-1; ++i){
    if (betas[i] <= t and t <= betas[i+1]){
      //std::cout << betas[i] << "   " << betas[i+1] << std::endl;
      return (rho[i+1]-rho[i])/(betas[i+1]-betas[i])*(t-betas[i])+rho[i];
    }
  }
  return 0.0;
}




template <typename Range>
void normalizeRho(const Range& betas, Range& rho){
  using namespace boost::math::quadrature;
  auto f = [betas,rho](const double& beta) { return getValue(betas,rho,beta); };
  double invArea = 1.0/(gauss<double,70>::integrate(f,0.0,betas[betas.size()-1]));
  for (auto& x : rho){ x *= invArea; }
}

template <typename Range>
auto getLambda_s(Range betas, Range rho){
  using std::sinh; using std::cosh; using std::exp;
  normalizeRho(betas, rho);

  std::vector<double> P (rho.size());
  for (size_t i = 1; i < rho.size(); ++i){
    P[i] = rho[i]/(2.0*betas[i]*sinh(betas[i]*0.5));
  }
  P[0] = rho[1]/(betas[1]*betas[1]);
  std::cout << std::endl;
  for ( auto x : rho){ std::cout << x << " ";}
  std::cout << std::endl;
  std::cout << std::endl;
  for ( auto x : P){ std::cout << x << " ";}
  std::cout << std::endl;

  using namespace boost::math::quadrature;
  auto integrand = [betas,rho](const double& beta) { 
    auto rhoVal = getValue(betas,rho,std::abs(beta));
    return (rhoVal/beta)*cosh(beta*0.5)/sinh(beta*0.5);
  };

  auto lambda_s = gauss<double,10>::integrate(integrand,0.0,betas[betas.size()-1]);

  return lambda_s;



}
