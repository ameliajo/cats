#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>

template <typename Range, typename Float>
auto sterp( Float B, Range betas, Range sLog ){
  Float SL;

  if      (B <= betas[0]             ){ return exp(sLog[0]); }
  else if (B >= betas[betas.size()-1]){ return 0.0;          }
  else {
    for ( size_t IC = 0; IC < betas.size()-1; ++IC ){
        if ( B >= betas[IC] and B <= betas[IC+1] ){
          SL = sLog[IC]+((B-betas[IC])/(betas[IC+1]-betas[IC]))*(sLog[IC+1]-sLog[IC]);
          return exp(SL);
        }
    }
  }
  return 0.0;
}





