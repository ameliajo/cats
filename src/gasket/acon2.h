#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>

template <typename Range, typename Float, typename RangeInts>
auto acon2(RangeInts NMAX, Range X5, Range ANK, Float T, Float SZCON, Range EPS, 
  Float AM, Float W1, Float PSQ, Range& S2) {
  using std::pow;

  int i, j, k, l, m, ind1, ind2, NE = S2.size();
  Range SK(1000);
  Float SINT, EIN, SINT_denominator = AM/(4.0*PSQ*T*W1), 
        PSQ_power = pow(PSQ*W1/AM,2);

  for ( size_t i = 0; i < S2.size(); ++i ){
    SK[i] = 0.0;
    l = NMAX[1]*2+1;
    for ( int j = 0; j < l; ++j ){
      m = NMAX[0]*2+1;
      ind2 = j-NMAX[1];
      for ( int k = 0; k < m; ++k){
        ind1 = k-NMAX[0];
        EIN = abs(EPS[i] - (k-NMAX[0])*X5[0] 
                         - (j-NMAX[1])*X5[1]);
        SINT = SZCON * exp(-(EIN*EIN+PSQ_power)*SINT_denominator);
        SK[i] = SK[i] + ANK[1+(abs(j-NMAX[1]))*X5.size()] * \
                        ANK[0+(abs(k-NMAX[0]))*X5.size()] * SINT;
      }
    }
    S2[i] = SK[i];
  }
}

