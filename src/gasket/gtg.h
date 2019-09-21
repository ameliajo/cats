#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>
#include "intg.h"

template <typename Range, typename Float, typename RangeInts>
auto GTG(Float wgt, Float T, Float AM, Range X, Range Q, Range t, Range PC, 
        Range PS, Float TBAR){
    int JS3 = X.size();
    int NT = t.size();
    int i;

//SUBROUTINE GTG(wgt,T,AM,X,Q,t,PC,PS,JS3,NT,TBAR)

//INTEGER :: JS3, NT
//REAL(8) :: X(JS3), Q(JS3), t(NT)
//REAL(8) :: X2(5), Q2(5), t2(10)
//REAL(8) :: wgt, T, AM, TBAR

    Float norm, F, H, C, CS, S;

    Float U = Q[0]*X[0]/3.0;
    // puts integral result in variable A -- JS3 is the length of X and Q
    //INTG(X,Q,A,JS3);
    double A = INTG(X,Q);

/*
norm = wgt / AM / (U+A)

CALL FTRANS(T,X,Q,t,PC,PS,JS3,NT)

F=X[0]*0.5/T
CALL COTH(H,F)
DO i=2,NT
  U=X[0]*t(i)
  if (U <= 0.005) {
    C=0.5*U - U**3/24.
    S=U-U**3/6.
    CS = U/3. - U**3/30.
  } else {
    C = COS(U)
    S = SIN(U)
    CS = S/U**2 - C/U
    C = (1.-C)/U
  }
  PC(i)=PC(i)+Q[0]/U*(H*(S-C)+C/F)
  PS(i)=PS(i)+Q[0]*CS
ENDDO
PC[0]= 0.5*Q[0]*(1./F+H)
PS[0]= 0.
TBAR = Q[0]*T*X[0]/3.

DO i=1,JS3
  F=X(i)*0.5/T
  CALL COTH(H,F)
  Q(i) = Q(i)*H/X(i)
ENDDO

CALL INTG(X,Q,A,JS3)
PC[0]=PC[0]+A

DO i=1,JS3
  Q(i) = Q(i)*X(i)**2*0.5
ENDDO
CALL INTG(X,Q,A,JS3)
TBAR= TBAR + A

DO i=1,NT
  PC(i) = PC(i)*norm
  PS(i) = PS(i)*norm
ENDDO

TBAR = TBAR*norm*AM/wgt

END SUBROUTINE
*/
}


