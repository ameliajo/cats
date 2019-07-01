#include <iostream>
#include <boost/math/quadrature/gauss.hpp>



template <typename Range, typename Float>
auto getValue(Range betas, Range rho, Float t){
  for (size_t i = 0; i < betas.size()-1; ++i){
    if (betas[i] <= t <= betas[i+1]){
      auto m = (rho[i+1]-rho[i])/(betas[i+1]-betas[i]);
      return m*(t-betas[i])+rho[i];
    }
  }
  return 0.0;
}



template <typename Range>
void normalizeRho(const Range& betas, Range& rho){
  double integral = 0.0;
  for (size_t i = 0; i < rho.size()-1; ++i){
    integral += (rho[i]+rho[i+1])*0.5*(betas[i+1]-betas[i]);
  }

  using namespace boost::math::quadrature;
  auto f = [betas,rho](const double& t) { return getValue(betas,rho,t); };
  double Q = gauss<double, 10>::integrate(f, 0.0, betas[betas.size()-1]);
  std::cout << integral<< std::endl;
  std::cout << Q<< std::endl;



  //double invArea = 1.0/integral;
  //for (auto& x : rho){
  //  x *= invArea;
  //}

}

auto quick(){
  using namespace boost::math::quadrature;
  auto f = [](const double& t) { return t * t * std::atan(t); };
  double Q = gauss<double, 7>::integrate(f, 0, 1);
  std::cout << Q<< std::endl;

}

template <typename Range>
auto getLambda_s(Range betas, Range rho){
  using std::sinh;
  using std::cosh;
  normalizeRho(betas, rho);
  return;

  auto c = rho[1]/(betas[1]*betas[1]);
  double integration = c*sinh(betas[1]*0.5);

  for (size_t i = 1; i < betas.size()-1; ++i){
    auto l = rho[i]/betas[i] * cosh(betas[i]*0.5) / sinh(betas[i]*0.5);
    auto r = rho[i+1]/betas[i+1] * cosh(betas[i+1]*0.5) / sinh(betas[i+1]*0.5);
    integration = integration + (betas[i+1]-betas[i])*0.5*(l+r);
  }

  quick();

  //std::cout << integration << std::endl;


}
