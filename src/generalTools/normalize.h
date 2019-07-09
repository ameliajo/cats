#ifndef CATS_NORMALIZE
#define CATS_NORMALIZE

#include <boost/math/quadrature/gauss.hpp>
#include "generalTools/interpolate.h"
template <typename Range>
void normalizeRho(const Range& betas, Range& rho, bool useOld){
  using namespace boost::math::quadrature;
  auto f1 = [betas,rho](auto i) { return rho[i]; };
  auto f2 = [betas,rho](const double& beta) { return interpolate(betas,rho,beta); };
  double invArea = useOld ? 
      1.0/(trapz(betas,f1))//-0.5*rho[1]*betas[1] + rho[1]*pow(betas[1],3)/3.0)
   :  1.0/(gauss<double,100>::integrate(f2,0.0,betas[betas.size()-1]));
  //std::cout << "area: "<< 1.0/invArea << std::endl;
  for (auto& x : rho){ x *= invArea; }
}

#endif
