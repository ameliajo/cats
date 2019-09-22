#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>
#include "sterp.h"

template <typename RangeInts, typename Float, typename Range>
auto rconv( RangeInts nMax, Range X5, Range ANK, Float T, Range& S1, Range betas){
  using std::abs; using std::log;

  int ind;
  Float SINT, BETAIN;
  Range SLOG(1000), SK(1000);

  for (size_t k = 0; k < X5.size(); ++k){
    for (size_t i = 0; i < S1.size(); ++i){
      if (S1[i] <= 0) { SLOG[i]=-100; } 
      else            { SLOG[i] = log(S1[i]); }
    }
    for ( size_t i = 0; i < S1.size(); ++i){
      SK[i]=0.0;
      for (int j = 0; j < nMax[k]*2+1; ++j){
        ind = j-nMax[k];
        BETAIN = abs(betas[i]-double(ind)*X5[k]/T);
        SINT = sterp(BETAIN,betas,SLOG);
        SK[i] += ANK[k+(abs(ind))*X5.size()] * SINT;
      }
    S1[i]=SK[i];
    }
  }
}


