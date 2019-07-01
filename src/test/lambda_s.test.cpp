#include "catch.hpp"
#include "getLambda_s.h"

TEST_CASE( "" ){
  GIVEN( "" ){ 
    REQUIRE( true );
    std::vector<double> rho   {0.0,1.0,2.0,3.0,3.5,5.0,6.0,8.0,6.0};
    std::vector<double> betas {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6};
    std::cout << std::endl;
    getLambda_s(betas,rho);
    std::cout << std::endl;

  } // GIVEN
} // TEST CASE
