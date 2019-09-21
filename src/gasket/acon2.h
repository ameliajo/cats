#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>

template <typename Range, typename Float, typename RangeInts>
auto acon2(RangeInts NMAX, Range X5, Range ANK, Float T, Float SZCON, Range EPS, 
  Float AM, Float W1, Float PSQ, Range S2) {

  int i, j, k, l, m, ind1, ind2, NE = S2.size();
  Range SK(1000);
  Float SINT, EIN;

  for ( size_t i = 0; i < S2.size(); ++i ){
    SK[i] = 0.0;
    for ( size_t j = 1; j < l; ++j ){
      m = NMAX[0]*2+1;
      for ( int k = 0; k < m; ){
        ind1 = k-NMAX[0]-1;
        EIN = abs(EPS[i]-ind1*X5[0]-ind2*X5[1]);
        SINT = SZCON*EXP(-AM*(EIN**2+(PSQ*W1/AM)**2)/(4.*PSQ*T*W1));
      }
    }
  }
  /*
  DO i=1,NE
    SK(i)=0.0
    l=NMAX(2)*2+1
    DO j=1,l
      m=NMAX(1)*2+1
      ind2 = j-NMAX(2)-1
      DO k=1,m
        ind1 = k-NMAX(1)-1
        EIN = ABS(EPS(i)-ind1*X5(1)-ind2*X5(2))
        SINT = SZCON*EXP(-AM*(EIN**2+(PSQ*W1/AM)**2)/(4.*PSQ*temperature*W1))
        if (SINT.GT.200) WRITE(*,*) SINT, ind1, ind2, ANK(2,ABS(ind2)+1), ANK(1,ABS(ind1)+1)
        SK(i) = SK(i) + ANK(2,ABS(ind2)+1)*ANK(1,ABS(ind1)+1)*SINT
      ENDDO
    ENDDO
    S2(i)=SK(i)
  ENDDO
END SUBROUTINE

*/
}

