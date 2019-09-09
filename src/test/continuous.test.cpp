#include "catch.hpp"
#include "continuous.h"

TEST_CASE( "" ){
  GIVEN( "simple H2O phonon distribution and its corresponding beta grid" ){ 
    double kb = 8.6173332e-5, T;
      std::vector<double> rho {0, .0005, .001, .002, .0035, .005, .0075, .01, 
        .013, .0165, .02, .0245, .029, .034, .0395, .045, .0506, .0562, .0622, 
      .0686, .075, .083, .091, .099, .107, .115, .1197, .1214, .1218, .1195, 
      .1125, .1065, .1005, .09542, .09126, .0871, .0839, .0807, .07798, .07574, 
      .0735, .07162, .06974, .06804, .06652, .065, .0634, .0618, .06022, .05866, 
      .0571, .05586, .05462, .0535, .0525, .0515, .05042, .04934, .04822, 
      .04706, .0459, .04478, .04366, .04288, .04244, .042, 0.};
    std::vector<double> betas (rho.size());
    T = 296.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.00255/(kb*T); }

    std::cout << getSab(betas,rho,0.01,0.006) << std::endl;
    std::cout << std::endl;
    //getSab(betas,rho,0.1,0.2);
  } // GIVEN
} // TEST CASE
/*
*/