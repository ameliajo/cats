#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>

template <typename RangeInts, typename Float, typename Range>
auto rconv( RangeInts nMax, Range X5, Range ANK, Float T, Range S1, Range betas){
    using std::abs; using std::log;
//SUBROUTINE RCONV(NE,JS5, NMAX,X5,ANK,temperature,S1,BETA)

//IMPLICIT NONE

    int i, j, k, ind, l;
    Float SINT, BETAIN;
    int NE = S1.size();
    int JS5 = X5.size();
    Range SLOG(1000), SK(1000);


//REAL(8) :: S1(NE), X5(JS5), temperature, BETA(NE),  ANK(JS5,20)
//REAL(8) :: SLOG(1000), SK(1000), SINT, BETAIN
  for (int k = 0; k < JS5; ++k){
    for (int i = 0; i < NE; ++i){
      if (S1[i] <= 0) {
        SLOG[i]=-100;
      } else {
        SLOG[i] = log(S1[i]);
      }
    }
    for ( int i = 0; i < NE; ++i){
      SK[i]=0.0;
      l=nMax[k]*2+1;
      for (int j = 0; j < l; ++j){
        ind = j-nMax[k]-1;
        BETAIN = abs(betas[i]-double(ind)*X5[k]/T);
        //std::cout << k << "    " << abs(ind)-1 << std::endl;
        //std::cout << ANK[k+(abs(ind)-1)*JS5] << std::endl;
        //std::cout << ANK[(k+3)+(abs(ind)-1+4)*JS5] << std::endl;
        //if (j == l-1 ) return;
  //    CALL STERP(BETAIN,BETA,NE,SINT,SLOG)
        SK[i] += ANK[k+(abs(ind)-1)*JS5] * SINT;
        }
    S1[i]=SK[i];
    }
  }
}


