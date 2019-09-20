#include "catch.hpp"
#include "../bessl.h"
#include "../../generalTools/testing.h"


TEST_CASE( "bessl" ){
  double xVal;
  std::vector<double> B (10,0.0), correctB (10);
  int nMax;

  GIVEN( "a small x value (xVal < 0.05)" ){
    WHEN( "results of one is fed into the next (round 1)" ){
      xVal = 7.4491602E-5;
      nMax = bessl(xVal,B,10);
      correctB = { 1, 3.72458E-5, 6.93624E-10, 8.61153E-15, 8.01859E-20, 0, 0, 0, 0, 0 };
      REQUIRE( ranges::equal(B,correctB,equal) );
      REQUIRE( nMax == 5 );

      xVal = 2.895559E-7;
      nMax = bessl(xVal,B,10);
      correctB = {1.0, 0.144778E-06, 0.104803E-13, 0, 0, 0, 0, 0, 0, 0};
      REQUIRE( ranges::equal(B,correctB,equal) );
      REQUIRE( nMax == 3 );
    } // WHEN
    WHEN( "results of one is fed into the next (round 2)" ){
      B = {1.0, 0.143330E-4, 0.102718E-9, 0.490751E-15, 0, 0, 0, 0, 0, 0};

      xVal = 7.4491602E-3;
      nMax = bessl(xVal,B,10);
      correctB = { 1.0, 0.372458E-2, 0.693625E-5, 0.861154E-8, 0.801859E-11,
                   0.597318E-14, 0.370793E-17, 0, 0, 0 };
      REQUIRE( ranges::equal(B,correctB,equal) );
      REQUIRE( nMax == 7 );

      xVal = 2.895559E-5;
      nMax = bessl(xVal,B,10);
      correctB = {1.0, 0.144778E-4, 0.104803E-9, 0.505773E-15, 0, 0, 0, 0, 0, 0};
      REQUIRE( ranges::equal(B,correctB,equal) );
      REQUIRE( nMax == 4 );
    } // WHEN
  } // GIVEN
  GIVEN( "a larger x value (xVal >= 0.05)" ){
      xVal = 0.05;
      nMax = bessl(xVal,B,10);
      correctB = {1.000625, 2.500781E-2, 3.125651E-4, 2.604574E-6, 1.627808E-8,
           8.138869E-11, 3.391145E-13, 1.211110E-15, 3.784685E-18, 1.051294E-20};
      REQUIRE( ranges::equal(B,correctB,equal) );
      REQUIRE( nMax == 0 );

      nMax = bessl(xVal,B,10);
      correctB = {1.000625, 2.500781E-2, 3.125651E-4, 2.604574E-6, 1.627808E-8,
      8.138869E-11, 3.391145E-13, 1.211110E-15, 3.784685E-18, 1.051294E-20};
      REQUIRE( ranges::equal(B,correctB,equal) );
      REQUIRE( nMax == 0 );


  } // GIVEN
} // TEST_CASE
      //std::cout << (B|ranges::view::all) << std::endl;
