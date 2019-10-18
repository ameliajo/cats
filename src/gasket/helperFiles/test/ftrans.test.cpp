#include "catch.hpp"
#include "../ftrans.h"
#include "../../../generalTools/testing.h"


TEST_CASE( "ftrans" ){
  GIVEN( "Simple, small X Q t vectors" ){
    std::vector<double> X(5), Q(5), correctPC(10), correctPS(10), t(10);
    double T = 2.55e-2;
    t = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
    {
      X = { 0.1, 0.2, 0.3, 0.5, 0.8 };
      Q = { 0.1, 0.3, 0.7, 0.4, 0.0 };
      auto out = ftrans(T, X, Q, t);
      correctPC = { 0.0, 0.752123, 0.752096, 0.752059, 0.752011, 0.751953, 
                    0.751883, 0.751804, 0.751713, 0.751612};
      correctPS = { 0.0, 0.00519993, 0.00779978, 0.0103995, 0.0129990, 
                    0.0155982, 0.0181972, 0.0207958, 0.0233940, 0.0259917};
      std::vector<double> PC = std::get<0>(out), PS = std::get<1>(out);

      REQUIRE( ranges::equal(PC,correctPC,equal) );
      REQUIRE( ranges::equal(PS,correctPS,equal) );
    }
    {
      X = { 0.1, 0.7, 0.10, 0.11, 0.20 };
      Q = { 0.0, 0.5, 0.0,  0.4,  0.0  };
      correctPC = { 0.0, 0.186750, 0.186749, 0.186748, 0.186746, 0.186745, 
                    0.186742, 0.186740, 0.186736, 0.186733 },
      correctPS = { 0.0, 0.496969E-3, 0.745452E-3, 0.993934E-3, 0.00124241,
                    0.00149089, 0.00173937, 0.00198784, 0.00223630, 0.00248477};
      auto out = ftrans(T, X, Q, t);
      std::vector<double> PC = std::get<0>(out), PS = std::get<1>(out);

      REQUIRE( ranges::equal(PC,correctPC,equal) );
      REQUIRE( ranges::equal(PS,correctPS,equal) );

    } // WHEN
  } // GIVEN
  GIVEN( "More realistic water example" ){
    std::vector<double> X(67), Q(67), correctPC(80), correctPS(80);
    Q = { 0.00049, 0.00098, 0.00191, 0.00338, 0.00485, 0.00721, 0.00966, 
          0.01245, 0.01588, 0.01931, 0.02353, 0.02794, 0.03260, 0.03799, 
          0.04338, 0.04884, 0.05433, 0.05994, 0.06622, 0.07249, 0.07971, 
          0.08755, 0.09539, 0.10324, 0.11108, 0.11730, 0.12050, 0.12158, 
          0.12081, 0.11662, 0.11015, 0.10426, 0.09838, 0.09403, 0.08995, 
          0.08616, 0.08302, 0.07988, 0.07745, 0.07526, 0.07313, 0.07129, 
          0.06945, 0.06783, 0.06634, 0.06484, 0.06327, 0.06171, 0.06016, 
          0.05863, 0.00000, 0.05588, 0.05467, 0.05356, 0.05258, 0.05160, 
          0.05055, 0.04949, 0.04840, 0.04726, 0.04613, 0.04502, 0.04392, 
          0.04299, 0.04256, 0.04213, 0.02471 };
    X = ranges::view::iota(1,int(Q.size()+1)) 
      | ranges::view::transform([](int i){return i*0.0025;});
    double T = 2.55e-2;
    std::vector<double> 
      t = ranges::view::iota(0,int(correctPC.size()))
        | ranges::view::transform([](int i){return i*0.1;});
    auto out = ftrans(T, X, Q, t);
    std::vector<double> PC = std::get<0>(out), PS = std::get<1>(out);

    correctPC = { 0.0, 0.194935, 0.194921, 0.194897, 0.194863, 0.194820, 
      0.194767, 0.194705, 0.194633, 0.194551, 0.194460, 0.194360, 0.194249, 
      0.194130, 0.194000, 0.193862, 0.193713, 0.193556, 0.193389, 0.193212, 
      0.193026, 0.192831, 0.192626, 0.192412, 0.192189, 0.191956, 0.191714, 
      0.191463, 0.191202, 0.190933, 0.190654, 0.190366, 0.190069, 0.189764, 
      0.189449, 0.189125, 0.188792, 0.188450, 0.188099, 0.187740, 0.187372, 
      0.186995, 0.186609, 0.186215, 0.185812, 0.185401, 0.184981, 0.184552, 
      0.184115, 0.183670, 0.183217, 0.182755, 0.182285, 0.181807, 0.181320, 
      0.180826, 0.180324, 0.179814, 0.179295, 0.178769, 0.178236, 0.177694, 
      0.177145, 0.176589, 0.176024, 0.175453, 0.174874, 0.174288, 0.173694, 
      0.173093, 0.172485, 0.171870, 0.171248, 0.170619, 0.169984, 0.169341, 
      0.168692, 0.168036, 0.167374, 0.166705 };
    correctPS = { 0.0, 0.977790E-3, 0.195549E-2, 0.293300E-2, 0.391024E-2,
      0.488711E-2, 0.586351E-2, 0.683937E-2, 0.781458E-2, 0.878906E-2,
      0.976270E-2, 0.0107354, 0.0117072, 0.0126778, 0.0136472, 0.0146153, 
      0.0155821, 0.0165474, 0.0175111, 0.0184732, 0.0194336, 0.0203921, 
      0.0213488, 0.0223034, 0.0232559, 0.0242063, 0.0251543, 0.0261001, 
      0.0270433, 0.0279841, 0.0289222, 0.0298577, 0.0307903, 0.0317201, 
      0.0326469, 0.0335706, 0.0344912, 0.0354087, 0.0363228, 0.0372335, 
      0.0381408, 0.0390445, 0.0399447, 0.0408411, 0.0417337, 0.0426224, 
      0.0435073, 0.0443880, 0.0452647, 0.0461372, 0.0470055, 0.0478694, 
      0.0487290, 0.0495840, 0.0504345, 0.0512803, 0.0521215, 0.0529579, 
      0.0537894, 0.0546160, 0.0554377, 0.0562542, 0.0570657, 0.0578719, 
      0.0586729, 0.0594686, 0.0602589, 0.0610437, 0.0618230, 0.0625968, 
      0.0633648, 0.0641272, 0.0648838, 0.0656346, 0.0663795, 0.0671184, 
      0.0678513, 0.0685782, 0.0692990, 0.0700136 };
    REQUIRE( ranges::equal(PC,correctPC,equal) );
    REQUIRE( ranges::equal(PS,correctPS,equal) );
  } // GIVEN

} // TEST_CASE