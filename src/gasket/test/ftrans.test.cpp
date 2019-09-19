#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../ftrans.h"

TEST_CASE( "ftrans" ){
  GIVEN( " " ){
    std::vector<double> 
        X { 0.1, 0.2, 0.3, 0.5, 0.8 }, 
        Q { 0.1, 0.3, 0.7, 0.4, 0.0 },
        t { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
    double T = 2.55e-2;
    WHEN( " " ){
        ftrans(T, X, Q, t);
        REQUIRE(true);
    } // WHEN
  } // GIVEN
} // TEST_CASE
