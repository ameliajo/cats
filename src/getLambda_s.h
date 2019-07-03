#include <iostream>
#include <boost/math/quadrature/gauss.hpp>

template <typename Range, typename Float>
auto getValue(Range betas, Range rho, Float t){
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
  auto f2 = [betas,rho](const double& beta) { return getValue(betas,rho,beta); };
  double invArea = useOld ? 
    1.0/trapz(betas,f1) : 
    1.0/(gauss<double,100>::integrate(f2,0.0,betas[betas.size()-1]));
  std::cout << "AN   " << 1.0/invArea << std::endl;
  for (auto& x : rho){ x *= invArea; }
}




template <typename Range>
auto getLambda_s(Range betas, Range rho, bool useOld = false){
  using std::sinh; using std::cosh; using std::exp;
  //normalizeRho(betas, rho, useOld);

  /*
  using namespace boost::math::quadrature;
  auto f1 = [betas,rho](int i) { 
    if (abs(betas[i]) < 1e-12){ 
      auto P = rho[1]/(betas[1]*betas[1]);
      return P*(exp(-betas[i]*0.5)+exp(betas[i]*0.5));
    }
    return (rho[i]/betas[i])*cosh(betas[i]*0.5)/sinh(betas[i]*0.5);
  };
  auto getP = [betas,rho](int i){ 
    if (i > 0){ return rho[i]/(2*betas[i]*sinh(betas[i]*0.5))*(exp(betas[i]*0.5)+exp(-betas[i]*0.5)); }
    return rho[1]/(betas[1]*betas[1]);
  };
  auto f2 = [betas,rho](const double& beta) { 
    if (beta < 1e-8){ 
      auto P = rho[1]/(betas[1]*betas[1]);
      return P*(exp(-beta*0.5)+exp(beta*0.5));
    }
    auto rhoVal = getValue(betas,rho,std::abs(beta));
    return (rhoVal/beta)*cosh(beta*0.5)/sinh(beta*0.5);
  };
  */
  double integral = 0.0;
  for (size_t i = 0; i < betas.size()-1; ++i){
    integral += 0.5*(rho[i]+rho[i+1])*(betas[i+1]-betas[i]);
  }
  integral = integral - 0.5*rho[1]*betas[1] + rho[1]*pow(betas[1],3)/3.0;

  for (auto& x : rho){ x /= integral; }


  std::vector<double> P(rho.size());
  P[0] = rho[1]/(betas[1]*betas[1]);
  for (size_t i = 1; i < P.size(); ++i){
    P[i] = rho[i]/(betas[i]*2*sinh(betas[i]*0.5));
  }

  std::cout << betas[0] << "    " << P[0] << std::endl;
  std::cout << betas[1] << "    " << P[1] << std::endl;
  std::cout << betas[2] << "    " << P[2] << std::endl;
  std::cout << betas[3] << "    " << P[3] << std::endl;

  normalizeRho(betas, P, useOld);


  //return useOld ? trapz(betas,f1)
  //              : gauss<double,10>::integrate(f2,0.0,betas[betas.size()-1]);
  return 1.0;

}
