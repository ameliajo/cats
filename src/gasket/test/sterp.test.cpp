#include "catch.hpp"
#include "../sterp.h"
#include "../../generalTools/testing.h"


TEST_CASE( "sterp" ){
  //double xVal;
  //std::vector<double> B (10,0.0), correctB (10);
  //int nMax;
    std::vector<double> betas {     0.000, 0.080, 0.160, 0.240, 0.320,
    0.400, 0.480, 0.560, 0.640, 0.720,
    0.800, 0.880, 0.960, 1.040, 1.120,
    1.200, 1.280, 1.360, 1.440, 1.520,
    1.600, 1.680, 1.760, 1.840, 1.920,
    2.000, 2.080, 2.160, 2.240, 2.320,
    2.400, 2.480, 2.560, 2.640, 2.720,
    2.800, 2.880, 2.960, 3.040, 3.120,
    3.200, 3.280, 3.360, 3.440, 3.520,
    3.600, 3.680, 3.760, 3.840, 3.920,
    4.000, 4.080, 4.160, 4.240, 4.320,
    4.400, 4.480, 4.560, 4.640, 4.720,
    4.800, 4.880, 4.960, 5.040, 5.120,
    5.200, 5.280, 5.360, 5.440, 5.520,
    5.600, 5.680, 5.760, 5.840, 5.920,
    6.000, 6.080, 6.160, 6.240, 6.320,
    6.400, 6.480, 6.560, 6.640, 6.720,
    6.800, 6.880, 6.960, 7.040, 7.120,
    7.200, 7.280, 7.360, 7.440, 7.520,
    7.600, 7.680, 7.760, 7.840, 7.920};

    double betaVal = 0.0;
    std::vector<double> nonZeroSLOGVals {     -2.8792469E+00, -2.9181933E+00, -2.9562901E+00, -2.9944738E+00, -3.0327443E+00,
    -3.0711013E+00, -3.1095448E+00, -3.1480745E+00, -3.1866902E+00, -3.2253920E+00,
    -3.2641795E+00, -3.3030527E+00, -3.3420115E+00, -3.3810558E+00, -3.4201855E+00,
    -3.4594004E+00, -3.4987006E+00, -3.5380859E+00, -3.5775563E+00, -3.6171117E+00,
    -3.6567521E+00, -3.6964774E+00, -3.7362876E+00, -3.7761828E+00, -3.8161628E+00,
    -3.8562277E+00, -3.8963775E+00, -3.9366122E+00, -3.9769318E+00, -4.0173364E+00,
    -4.0578260E+00, -4.0984007E+00, -4.1390604E+00, -4.1798053E+00, -4.2206355E+00,
    -4.2615510E+00, -4.3025520E+00, -4.3436385E+00, -4.3848106E+00, -4.4260684E+00,
    -4.4674122E+00, -4.5088419E+00, -4.5503579E+00, -4.5919601E+00, -4.6336489E+00,
    -4.6754243E+00, -4.7172865E+00, -4.7592358E+00, -4.8012723E+00, -4.8433962E+00,
    -4.8856079E+00, -4.9279074E+00, -4.9702951E+00, -5.0127711E+00, -5.0553358E+00,
    -5.0979895E+00, -5.1407323E+00, -5.1835646E+00, -5.2264868E+00, -5.2694991E+00,
    -5.3126018E+00, -5.3557953E+00, -5.3990799E+00, -5.4424560E+00, -5.4859240E+00,
    -5.5294843E+00, -5.5731371E+00, -5.6168831E+00, -5.6607225E+00, -5.7046558E+00,
    -5.7486836E+00, -5.7928061E+00, -5.8370239E+00, -5.8813375E+00, -5.9257474E+00,
    -5.9702542E+00, -6.0148583E+00, -6.0595602E+00, -6.1043607E+00, -6.1492602E+00,
    -6.1942593E+00, -6.2393587E+00, -6.2845590E+00, -6.3298608E+00, -6.3752648E+00,
    -6.4207717E+00, -6.4663822E+00, -6.5120971E+00, -6.5579170E+00, -6.6038427E+00,
    -6.6498750E+00, -6.6960148E+00, -6.7422628E+00, -6.7886199E+00, -6.8350870E+00,
    -6.8816649E+00, -6.9283547E+00, -6.9751571E+00, -7.0220734E+00, -7.0691043E+00};
    std::vector<double> sLog (1000,0.0);
    for (size_t i = 0; i < nonZeroSLOGVals.size(); ++i){
        sLog[i] = nonZeroSLOGVals[i];
    }
    
    REQUIRE( 5.61770521E-2 == Approx(sterp(0.00, betas, sLog)).epsilon(1e-6) );
    REQUIRE( 5.40312183E-2 == Approx(sterp(0.08, betas, sLog)).epsilon(1e-6) );
    REQUIRE( 5.20115180E-2 == Approx(sterp(0.16, betas, sLog)).epsilon(1e-6) );
    REQUIRE( 5.17638601E-2 == Approx(sterp(0.17, betas, sLog)).epsilon(1e-6) );
    REQUIRE( 8.51495284E-4 == Approx(sterp(7.919, betas, sLog)).epsilon(1e-6) );
    REQUIRE( 0.0 == Approx(sterp(7.920, betas, sLog)).epsilon(1e-6) );
    REQUIRE( 5.61770521E-2 == Approx(sterp(-0.01, betas, sLog)).epsilon(1e-6) );

  GIVEN( "" ){
    WHEN( "" ){

    } // WHEN
  } // GIVEN
} // TEST_CASE
      //std::cout << (B|ranges::view::all) << std::endl;
