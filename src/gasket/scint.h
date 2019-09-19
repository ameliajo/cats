


template <typename Range, typename Float>
auto scint(Range t,Range GC,Range GS,Range EPS,Range S, Float T, Float A, Float B , Float F, int NT, int NE){
    /*
    INTEGER :: NT, NE
    REAL(8) :: t(NT), GC(NT), GS(NT), S(NE),EPS(NE)
    REAL(8) :: A, B, F, temperature

    INTEGER :: i, j

    REAL(8) :: SM, SINSM, COSSM, V0, U, SINS, COSS, V, ST, CT, SINT, COST, AL
    REAL(8) :: EX1, EX2
    REAL(8) :: Q1,Q2,R1,R2
    */

EX1 = EXP(B*GC[0])
EX2 = EXP(-A*temperature*t[0]*t[0])
Q1 = (COS(B*GS[0])*EX1-1.)*EX2
R1 = SIN(B*GS[0])*EX1*EX2

DO i=1,NE
  AL = A-EPS[i]
  IF (AL.EQ.0) THEN
    WRITE(*,*) 'AL IS EQUAL to 0'
  ELSE
    S(i) = 0.
    SM = t[0]*AL
    SINSM = SIN(SM)
    COSSM = COS(SM)
    V0 = 0.
    DO j=2,NT
      U=t[j]*AL
      SINS = SIN(U)
      COSS = COS(U)
      V=U-SM
      IF (ABS(V/V0-1.).LE.5e-7) GO TO 40
      IF (ABS(V).LE.0.005) THEN
        ST = (V*V)/6.-(V*V*V*V)/120.
        CT = V*0.5-V*V*V/24.
      ELSE
        SINT = SINS*COSSM-COSS*SINSM
        COST = COSS*COSSM+SINS*SINSM
        ST = 1.-SINT/V
        CT = (1.-COST)/V
      ENDIF
  40  EX1 = EXP(B*GC(j))
      EX2 = EXP(-A*temperature*t[j]*t[j])
      Q2 = (COS(B*GS(j))*EX1-1.)*EX2
      R2 = SIN(B*GS(j))*EX1*EX2
      S(i)=S(i)+Q2*(ST*SINS+CT*COSS)-Q1*(ST*SINSM-CT*COSSM)
      S(i)=S(i)-(R2*(CT*SINS-ST*COSS)+R1*(ST*COSSM+CT*SINSM))
      SM=U
      SINSM=SINS
      COSSM=COSS
      V0=V
      Q1=Q2
      R1=R2
    ENDDO
    S(i)=S(i)*F/AL
  ENDIF
ENDDO

