#include "catch.hpp"
#include "getLambda_s.h"

TEST_CASE( "Debye waller factor, lambda_s (trapezoid integral)" ){
  double kb = 8.6173332e-5, T;
  GIVEN( "simple H2O phonon distribution and its corresponding beta grid" ){ 
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
    auto lambda_s = getLambda_s(betas,rho,true);
    REQUIRE( 0.52920997122117586 == Approx(lambda_s).epsilon(1e-6) );

    T = 600.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.00255/(kb*T); }
    lambda_s = getLambda_s(betas,rho,true);
    REQUIRE( 1.760282876 == Approx(lambda_s).epsilon(1e-6) );

  } // GIVEN

  GIVEN( "Be in BeO phonon distribution and its corresponding beta grid" ){ 
    std::vector<double> rho { 3.100000E-6, 6.899985E-6, 1.579987E-5, 3.249975E-5, 
    5.624949E-5, 8.334943E-5, 1.141990E-4, 1.542986E-4, 1.998983E-4, 2.561974E-4, 
    3.232470E-4, 3.888969E-4, 4.642953E-4, 5.549946E-4, 6.520937E-4, 7.426447E-4, 
    8.467409E-4, 9.801895E-4, 1.118239E-3, 1.260638E-3, 1.427084E-3, 1.598535E-3, 
    1.781880E-3, 1.997277E-3, 2.222176E-3, 2.447475E-3, 2.693269E-3, 2.959668E-3, 
    3.249061E-3, 3.548913E-3, 3.871502E-3, 4.231850E-3, 4.647983E-3, 5.072545E-3, 
    5.508525E-3, 5.996425E-3, 6.571393E-3, 7.270280E-3, 8.078150E-3, 9.242100E-3, 
    1.058603E-2, 1.170882E-2, 1.260280E-2, 1.355699E-2, 1.489420E-2, 1.686302E-2, 
    1.985052E-2, 2.108993E-2, 1.985079E-2, 1.877167E-2, 1.770130E-2, 1.647055E-2, 
    1.539779E-2, 1.454566E-2, 1.401839E-2, 1.371271E-2, 1.347211E-2, 1.325255E-2, 
    1.314101E-2, 1.313275E-2, 1.303920E-2, 1.263506E-2, 1.216104E-2, 1.174264E-2, 
    1.108204E-2, 1.048161E-2, 1.010476E-2, 9.748349E-3, 9.503098E-3, 9.594495E-3, 
    9.700141E-3, 9.128369E-3, 8.656028E-3, 9.131014E-3, 9.605228E-3, 9.105205E-3, 
    8.259261E-3, 8.786475E-3, 1.086055E-2, 1.277622E-2, 1.621424E-2, 1.946884E-2, 
    2.976416E-2, 9.697986E-2, 1.552680E-1, 1.328978E-1, 1.235346E-1, 1.214221E-1, 
    8.868644E-2, 5.501578E-2, 3.047213E-2, 1.783501E-2, 8.844040E-3, 4.163439E-3, 
    2.958511E-3, 2.121391E-3, 2.761060E-3, 5.546849E-3, 9.687513E-3, 1.308667E-2, 
    1.468507E-2, 1.572126E-2, 1.673388E-2, 1.808334E-2, 2.032280E-2, 2.189195E-2, 
    2.199840E-2, 2.159049E-2, 2.082654E-2, 2.048958E-2, 2.079742E-2, 2.099254E-2, 
    2.091494E-2, 1.985058E-2, 1.765808E-2, 1.589056E-2, 1.491489E-2, 1.428881E-2, 
    1.395833E-2, 1.452532E-2, 1.367085E-2, 1.145167E-2, 1.047933E-2, 1.045589E-2, 
    1.133254E-2, 1.124611E-2, 9.129051E-3, 7.375213E-3, 6.719129E-3, 6.179936E-3, 
    5.740469E-3, 5.386990E-3, 5.107686E-3, 4.892018E-3, 4.627346E-3, 4.353081E-3, 
    4.173784E-3, 4.024547E-3, 3.690213E-3, 3.332776E-3, 3.125731E-3, 2.892161E-3, 
    2.661229E-3, 2.460275E-3, 2.270017E-3, 1.943553E-3, 1.456127E-3, 9.878817E-4, 
    3.999894E-4, 1.361800E-5 };

    std::vector<double> betas (rho.size());

    T = 296.3;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 0.457453009 == Approx(getLambda_s(betas,rho,true)).epsilon(1e-6) );

    T = 600.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 1.451808074 == Approx(getLambda_s(betas,rho,true)).epsilon(1e-6) );

    T = 900.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 3.07133422 == Approx(getLambda_s(betas,rho,true)).epsilon(1e-6) );

  } // GIVEN
  GIVEN( "Crystalline graphite phonon distribution w/ corresponding beta grid" ){ 
    std::vector<double> rho { 0.000000E+0,2.389404E-4,6.138146E-4,1.092778E-3,
      1.608394E-3,2.186662E-3,2.803740E-3,3.272769E-3,3.823771E-3,4.337184E-3,
      4.965104E-3,5.504338E-3,6.226192E-3,6.904866E-3,7.907053E-3,8.928416E-3,
      1.009259E-2,1.175593E-2,1.385391E-2,1.692715E-2,2.110828E-2,1.684029E-2,
      1.544876E-2,1.469588E-2,1.423131E-2,1.405989E-2,1.377400E-2,1.358826E-2,
      1.389902E-2,1.372978E-2,1.366908E-2,1.397781E-2,1.401767E-2,1.413162E-2,
      1.430690E-2,1.464399E-2,1.489201E-2,1.513803E-2,1.515254E-2,1.562486E-2,
      1.583616E-2,1.605693E-2,1.655679E-2,1.662095E-2,1.698576E-2,1.763674E-2,
      1.748583E-2,1.823559E-2,1.867873E-2,1.907488E-2,1.970051E-2,2.047309E-2,
      2.104205E-2,2.178073E-2,2.270122E-2,2.408677E-2,2.603356E-2,2.856833E-2,
      3.525877E-2,3.908980E-2,3.102058E-2,1.783513E-2,1.566200E-2,1.462207E-2,
      1.344832E-2,1.263522E-2,1.177404E-2,1.114559E-2,1.178113E-2,1.345275E-2,
      1.530459E-2,1.747794E-2,1.953481E-2,2.207491E-2,2.520249E-2,2.988527E-2,
      3.742777E-2,4.957316E-2,5.206670E-2,4.752219E-2,4.049131E-2,3.708749E-2,
      3.515284E-2,3.309387E-2,3.262886E-2,3.194211E-2,3.097576E-2,3.058658E-2,
      3.035237E-2,2.979172E-2,2.967801E-2,2.933133E-2,2.951316E-2,2.909084E-2,
      2.926866E-2,2.911070E-2,2.878851E-2,2.855165E-2,2.953627E-2,2.959317E-2,
      2.923019E-2,2.981201E-2,3.090659E-2,3.137655E-2,3.305021E-2,3.454927E-2,
      3.662437E-2,4.033480E-2,2.179915E-2,1.058875E-2,1.032820E-2,1.007510E-2,
      1.004608E-2,1.011297E-2,9.763660E-3,9.915683E-3,9.485471E-3,1.009081E-2,
      1.001289E-2,9.846066E-3,9.893285E-3,9.769748E-3,9.878086E-3,9.692190E-3,
      6.749098E-3,7.110376E-3,7.415005E-3,7.517847E-3,7.538432E-3,7.734181E-3,
      7.888265E-3,8.211207E-3,8.189280E-3,8.591909E-3,8.780248E-3,8.889933E-3,
      9.266771E-3,9.277766E-3,9.477660E-3,1.010152E-2,1.010663E-2,1.040913E-2,
      1.101126E-2,1.114553E-2,1.145906E-2,1.177827E-2,1.229193E-2,1.287069E-2,
      1.334672E-2,1.432119E-2,1.584815E-2,1.461848E-2,1.164786E-2,1.150899E-2,
      1.147745E-2,1.179567E-2,1.195104E-2,1.168973E-2,1.198932E-2,1.212189E-2,
      1.212663E-2,1.202596E-2,1.245059E-2,1.222553E-2,1.266497E-2,1.194287E-2,
      1.104797E-2,1.328359E-2,1.457517E-2,4.796585E-2,6.635506E-2,5.952822E-2,
      5.822476E-2,6.230171E-2,7.723869E-2,7.414227E-2,4.922440E-2,4.101475E-2,
      3.677530E-2,3.381587E-2,3.225044E-2,3.096246E-2,2.980026E-2,2.865426E-2,
      2.839644E-2,2.789940E-2,2.775656E-2,2.792887E-2,2.774652E-2,2.790451E-2,
      2.826953E-2,2.853382E-2,2.938374E-2,3.067675E-2,3.248671E-2,3.480780E-2,
      3.904353E-2,5.496510E-2,7.153700E-2,2.967863E-2};

    std::vector<double> betas (rho.size());

    T = 296.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 0.8680557717 == Approx(getLambda_s(betas,rho,true)).epsilon(5e-3) );

    T = 600.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 3.1892824506 == Approx(getLambda_s(betas,rho,true)).epsilon(5e-3) );

    T = 900.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 6.9917441072 == Approx(getLambda_s(betas,rho,true)).epsilon(5e-3) );

  } // GIVEN

} // TEST CASE

TEST_CASE( "Debye waller factor, lambda_s (quadrature)" ){
  double kb = 8.6173332e-5, T;
  GIVEN( "simple H2O phonon distribution and its corresponding beta grid" ){ 
    std::vector<double> rho {0, .0005, .001, .002, .0035, .005, .0075, .01, 
      .013, .0165, .02, .0245, .029, .034, .0395, .045, .0506, .0562, .0622, 
      .0686, .075, .083, .091, .099, .107, .115, .1197, .1214, .1218, .1195, 
      .1125, .1065, .1005, .09542, .09126, .0871, .0839, .0807, .07798, .07574, 
      .0735, .07162, .06974, .06804, .06652, .065, .0634, .0618, .06022, .05866, 
      .0571, .05586, .05462, .0535, .0525, .0515, .05042, .04934, .04822, 
      .04706, .0459, .04478, .04366, .04288, .04244, .042, 0.};
    std::vector<double> betas (rho.size());

    /*
    T = 296.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.00255/(kb*T); }
    auto lambda_s = getLambda_s(betas,rho);
    REQUIRE( 0.529209971 == Approx(lambda_s).epsilon(1e-1) );

    T = 600.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.00255/(kb*T); }
    lambda_s = getLambda_s(betas,rho);
    REQUIRE( 1.760282876 == Approx(lambda_s).epsilon(1e-1) );
    */

  } // GIVEN
  GIVEN( "Be in BeO phonon distribution and its corresponding beta grid" ){ 
    std::vector<double> rho { 3.100000E-6, 6.899985E-6, 1.579987E-5, 3.249975E-5, 
    5.624949E-5, 8.334943E-5, 1.141990E-4, 1.542986E-4, 1.998983E-4, 2.561974E-4, 
    3.232470E-4, 3.888969E-4, 4.642953E-4, 5.549946E-4, 6.520937E-4, 7.426447E-4, 
    8.467409E-4, 9.801895E-4, 1.118239E-3, 1.260638E-3, 1.427084E-3, 1.598535E-3, 
    1.781880E-3, 1.997277E-3, 2.222176E-3, 2.447475E-3, 2.693269E-3, 2.959668E-3, 
    3.249061E-3, 3.548913E-3, 3.871502E-3, 4.231850E-3, 4.647983E-3, 5.072545E-3, 
    5.508525E-3, 5.996425E-3, 6.571393E-3, 7.270280E-3, 8.078150E-3, 9.242100E-3, 
    1.058603E-2, 1.170882E-2, 1.260280E-2, 1.355699E-2, 1.489420E-2, 1.686302E-2, 
    1.985052E-2, 2.108993E-2, 1.985079E-2, 1.877167E-2, 1.770130E-2, 1.647055E-2, 
    1.539779E-2, 1.454566E-2, 1.401839E-2, 1.371271E-2, 1.347211E-2, 1.325255E-2, 
    1.314101E-2, 1.313275E-2, 1.303920E-2, 1.263506E-2, 1.216104E-2, 1.174264E-2, 
    1.108204E-2, 1.048161E-2, 1.010476E-2, 9.748349E-3, 9.503098E-3, 9.594495E-3, 
    9.700141E-3, 9.128369E-3, 8.656028E-3, 9.131014E-3, 9.605228E-3, 9.105205E-3, 
    8.259261E-3, 8.786475E-3, 1.086055E-2, 1.277622E-2, 1.621424E-2, 1.946884E-2, 
    2.976416E-2, 9.697986E-2, 1.552680E-1, 1.328978E-1, 1.235346E-1, 1.214221E-1, 
    8.868644E-2, 5.501578E-2, 3.047213E-2, 1.783501E-2, 8.844040E-3, 4.163439E-3, 
    2.958511E-3, 2.121391E-3, 2.761060E-3, 5.546849E-3, 9.687513E-3, 1.308667E-2, 
    1.468507E-2, 1.572126E-2, 1.673388E-2, 1.808334E-2, 2.032280E-2, 2.189195E-2, 
    2.199840E-2, 2.159049E-2, 2.082654E-2, 2.048958E-2, 2.079742E-2, 2.099254E-2, 
    2.091494E-2, 1.985058E-2, 1.765808E-2, 1.589056E-2, 1.491489E-2, 1.428881E-2, 
    1.395833E-2, 1.452532E-2, 1.367085E-2, 1.145167E-2, 1.047933E-2, 1.045589E-2, 
    1.133254E-2, 1.124611E-2, 9.129051E-3, 7.375213E-3, 6.719129E-3, 6.179936E-3, 
    5.740469E-3, 5.386990E-3, 5.107686E-3, 4.892018E-3, 4.627346E-3, 4.353081E-3, 
    4.173784E-3, 4.024547E-3, 3.690213E-3, 3.332776E-3, 3.125731E-3, 2.892161E-3, 
    2.661229E-3, 2.460275E-3, 2.270017E-3, 1.943553E-3, 1.456127E-3, 9.878817E-4, 
    3.999894E-4, 1.361800E-5 };

    std::vector<double> betas (rho.size());

    /*
    T = 296.3;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 0.457453009 == Approx(getLambda_s(betas,rho,true)).epsilon(1e-6) );

    T = 600.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 1.451808074 == Approx(getLambda_s(betas,rho,true)).epsilon(1e-6) );

    T = 900.0;
    for (size_t i = 0; i < betas.size(); ++i){ betas[i] = i*0.001/(kb*T); }
    REQUIRE( 3.07133422 == Approx(getLambda_s(betas,rho,true)).epsilon(1e-5) );


    */

  } // GIVEN
} // TEST CASE
