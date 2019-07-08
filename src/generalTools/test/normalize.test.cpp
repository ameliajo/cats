#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "getLambda_s.h"

TEST_CASE( "Normalizing rho" ){
  using std::pow;
  using std::cos;

  std::vector<double> rho (1001);
  std::vector<double> betas (1001);

  GIVEN( "Parabolic vector spanning 0->1" ){ 
    auto parabola = [](const auto& beta){ return -beta*beta + beta; };
    for (size_t i = 0; i < rho.size(); ++i){ 
      betas[i] = double(i)/(rho.size()-1);
      rho[i]   = parabola(betas[i]);
    }
    THEN( "Result nearly equals analytic solution" ){
      double sum = 0.1666666667;
      normalizeRho(betas,rho,false);
      for (size_t i = 0; i < rho.size(); ++i){
        REQUIRE( rho[i] == Approx(parabola(betas[i])/sum).epsilon(1e-5) );
      }
    } // THEN
  } // GIVEN

  GIVEN( "Fourth order polynomial vector spanning 0->2" ){ 
    auto f = [](const auto beta){ return -pow(beta,4) + pow(beta,3); };
    for (size_t i = 0; i < rho.size(); ++i){ 
      betas[i] = 2.0*double(i)/(rho.size()-1);
      rho[i]   = f(betas[i]);
    }
    THEN( "Result nearly equals analytic solution" ){
      double sum = -2.4;
      normalizeRho(betas,rho,false);
      for (size_t i = 0; i < rho.size(); ++i){
        REQUIRE( rho[i] == Approx(f(betas[i])/sum).epsilon(1e-5) );
      }
    } // THEN
  } // GIVEN

  GIVEN( "Larger equation spanning 0->10" ){ 
    auto f = [](const auto beta){return -pow(beta,1.5) + cos(beta*3.0)*beta; };
    for (size_t i = 0; i < rho.size(); ++i){ 
      betas[i] = 10.0*double(i)/(rho.size()-1);
      rho[i]   = f(betas[i]);
    }
    THEN( "Result nearly equals analytic solution" ){
      double sum = -129.879;
      normalizeRho(betas,rho,false);
      for (size_t i = 0; i < rho.size(); ++i){
        REQUIRE( rho[i] == Approx(f(betas[i])/sum).epsilon(1e-5) );
      }
    } // THEN
  } // GIVEN
} // TEST CASE
