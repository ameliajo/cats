PROGRAM MAIN

IMPLICIT NONE

INTEGER :: nTime, nBeta
INTEGER :: nRho, JS4, nOsc
REAL(8) :: W1, W2, W3, W4, W5, temp

INTEGER :: i, j, ialpha

REAL(8) :: AM, EX, SINHEX, TBAR, SZCON

REAL(8), ALLOCATABLE :: X(:), Q(:), ZPHON(:), T(:), GC(:), GS(:)

REAL(8) :: ALPHA, PSQ !allocatable eventually
REAL(8) :: DBW, DBWP, APS, BPS, EX1, EX2, Q5(2), X5(2), ARG1(2), ARG2(2)
REAL(8), ALLOCATABLE :: S1(:), S2(:), S(:), BETAS(:), EPS(:), ALPHAS(:)
REAL(8) :: RR, ANK(2,10), U, BF(10), tstep
INTEGER :: NMAX(2), nAlpha

nTime = 700  !number of time steps
tstep = 0.1  !time step size

nBeta = 10 
nAlpha = 20

JS4 = 0
nOsc = 2   !number of discrete osciallators

Q5(1) = 0.33333   !discrete oscillators weights
Q5(2) = 0.66667

X5(1) = 0.205     !discrete oscillators energy
X5(2) = 0.480

W1 = 0.0555556    !free gas contribution
W2 = 0.
W3 = 0.4444444    !continuous spectrum
W4 = 0.
W5 = 0.5          !discrete osciallators

!Temperatures in eV
temp = 0.0255

OPEN(10, FILE='out2.txt')

!atomic mass
AM = 1.0/1.0086654 !Convert mass to neutron mass unit

nRho = 67
ALLOCATE(X(nRho), Q(nRho), ZPHON(nRho), ALPHAS(nAlpha))
X(1) = 0.0025
X(2) = 0.005
X(3) = 0.0075
X(4) = 0.01
X(5) = 0.0125
X(6) = 0.015
X(7) = 0.0175
X(8) = 0.02
X(9) = 0.0225
X(10) = 0.025
X(11) = 0.0275
X(12) = 0.03
X(13) = 0.0325
X(14) = 0.035
X(15) = 0.0375
X(16) = 0.04
X(17) = 0.0425
X(18) = 0.045
X(19) = 0.0475
X(20) = 0.05
X(21) = 0.0525
X(22) = 0.055
X(23) = 0.0575
X(24) = 0.06
X(25) = 0.0625
X(26) = 0.065
X(27) = 0.0675
X(28) = 0.07
X(29) = 0.0725
X(30) = 0.075
X(31) = 0.0775
X(32) = 0.08
X(33) = 0.0825
X(34) = 0.085
X(35) = 0.0875
X(36) = 0.09
X(37) = 0.0925
X(38) = 0.095
X(39) = 0.0975
X(40) = 0.1
X(41) = 0.1025
X(42) = 0.105
X(43) = 0.1075
X(44) = 0.11
X(45) = 0.1125
X(46) = 0.115
X(47) = 0.1175
X(48) = 0.12
X(49) = 0.1225
X(50) = 0.125
X(51) = 0.1275
X(52) = 0.13
X(53) = 0.1325
X(54) = 0.135
X(55) = 0.1375
X(56) = 0.14
X(57) = 0.1425
X(58) = 0.145
X(59) = 0.1475
X(60) = 0.15
X(61) = 0.1525
X(62) = 0.155
X(63) = 0.1575
X(64) = 0.16
X(65) = 0.1625
X(66) = 0.165
X(67) = 0.1675
!X(1) = 0.003;  X(2) = 0.007;  X(3) = 0.01;   X(4) = 0.014;  X(5) = 0.017; 
!X(6) = 0.021;  X(7) = 0.024;  X(8) = 0.028;  X(9) = 0.031;  X(10) = 0.035; 
!X(11) = 0.038; X(12) = 0.042; X(13) = 0.045; X(14) = 0.049; X(15) = 0.052; 
!X(16) = 0.056; X(17) = 0.059; X(18) = 0.062; X(19) = 0.066; X(20) = 0.069; 
!X(21) = 0.073; X(22) = 0.076; X(23) = 0.08;  X(24) = 0.083; X(25) = 0.087; 
!X(26) = 0.09;  X(27) = 0.094; X(28) = 0.097; X(29) = 0.101; X(30) = 0.104; 
!X(31) = 0.108; X(32) = 0.111; X(33) = 0.114; X(34) = 0.118; X(35) = 0.121; 
!X(36) = 0.125; X(37) = 0.128; X(38) = 0.132; X(39) = 0.135; X(40) = 0.139; 
!X(41) = 0.142; X(42) = 0.146; X(43) = 0.149; X(44) = 0.153; X(45) = 0.156; 
!X(46) = 0.16;  X(47) = 0.163; X(48) = 0.167; 

!Q(1)  = 0.001; Q(2)  = 0.002; Q(3)  = 0.004; Q(4) = 0.006;  Q(5)  = 0.01; 
!Q(6)  = 0.014; Q(7)  = 0.018; Q(8)  = 0.024; Q(9) = 0.03;   Q(10) = 0.037; 
!Q(11) = 0.045; Q(12) = 0.052; Q(13) = 0.06;  Q(14) = 0.069; Q(15) = 0.078; 
!Q(16) = 0.089; Q(17) = 0.1;   Q(18) = 0.111; Q(19) = 0.119; Q(20) = 0.121; 
!Q(21) = 0.12;  Q(22) = 0.113; Q(23) = 0.105; Q(24) = 0.097; Q(25) = 0.091; 
!Q(26) = 0.086; Q(27) = 0.082; Q(28) = 0.078; Q(29) = 0.075; Q(30) = 0.072; 
!Q(31) = 0.069; Q(32) = 0.067; Q(33) = 0.065; Q(34) = 0.063; Q(35) = 0.061; 
!Q(36) = 0.059; Q(37) = 0.057; Q(38) = 0.055; Q(39) = 0.053; Q(40) = 0.052; 
!Q(41) = 0.051; Q(42) = 0.049; Q(43) = 0.048; Q(44) = 0.046; Q(45) = 0.045; 
!Q(46) = 0.043; Q(47) = 0.042; Q(48) = 0.034; 
Q(1) = 0.00049
Q(2) = 0.00098
Q(3) = 0.00191
Q(4) = 0.00338
Q(5) = 0.00485
Q(6) = 0.00721
Q(7) = 0.00966
Q(8) = 0.01245
Q(9) = 0.01588
Q(10) = 0.01931
Q(11) = 0.02353
Q(12) = 0.02794
Q(13) = 0.0326
Q(14) = 0.03799
Q(15) = 0.04338
Q(16) = 0.04884
Q(17) = 0.05433
Q(18) = 0.05994
Q(19) = 0.06622
Q(20) = 0.07249
Q(21) = 0.07971
Q(22) = 0.08755
Q(23) = 0.09539
Q(24) = 0.10324
Q(25) = 0.11108
Q(26) = 0.1173
Q(27) = 0.1205
Q(28) = 0.12158
Q(29) = 0.12081
Q(30) = 0.11662
Q(31) = 0.11015
Q(32) = 0.10426
Q(33) = 0.09838
Q(34) = 0.09403
Q(35) = 0.08995
Q(36) = 0.08616
Q(37) = 0.08302
Q(38) = 0.07988
Q(39) = 0.07745
Q(40) = 0.07526
Q(41) = 0.07313
Q(42) = 0.07129
Q(43) = 0.06945
Q(44) = 0.06783
Q(45) = 0.06634
Q(46) = 0.06484
Q(47) = 0.06327
Q(48) = 0.06171
Q(49) = 0.06016
Q(50) = 0.05863
Q(51) = 0.0
Q(52) = 0.05588
Q(53) = 0.05467
Q(54) = 0.05356
Q(55) = 0.05258
Q(56) = 0.0516
Q(57) = 0.05055
Q(58) = 0.04949
Q(59) = 0.0484
Q(60) = 0.04726
Q(61) = 0.04613
Q(62) = 0.04502
Q(63) = 0.04392
Q(64) = 0.04299
Q(65) = 0.04256
Q(66) = 0.04213
Q(67) = 0.02471

DO i=1,nRho
  EX=EXP(X(I)/(2.*temp))
  SINHEX=0.5*(EX-1./EX)      ! not needed
  ZPHON(I) = Q(I)*EX/(2.*X(I)*SINHEX)  !not needed
ENDDO

!setup time grid
ALLOCATE(t(nTime))
t(1) = 0.0
DO i=2,nTime
  !t(i) = t(i-1)+tstep
  IF (i.GT.200.AND.i.LE.400) THEN
     t(i) = t(i-1)+tstep*5
  ELSEIF (i.GT.400.AND.i.LE.600) THEN
     t(i) = t(i-1)+tstep*10
  ELSEIF (i.GT.600.AND.i.LE.800) THEN
     t(i) = t(i-1)+tstep*20
  ELSEIF (i.GT.800.AND.i.LE.1000) THEN
     t(i) = t(i-1)+tstep*50
  ELSE
     t(i) = t(i-1)+tstep
  END IF
ENDDO


ALLOCATE(GC(nTime),GS(nTime))

CALL GTG(W3,temp,AM,X,Q,t,GC,GS,nRho,nTime,TBAR)

ALLOCATE(S1(nBeta),S2(nBeta),S(nBeta))

!setup beta grid
ALLOCATE(BETAS(nBeta),EPS(nBeta))
BETAS(1) = 0.0
EPS(1) = 0.0
DO i=2,nBeta
  BETAS(i) = BETAS(i-1)+0.08
  EPS(i) = BETAS(i)*temp
ENDDO

ALPHAS(1) = 0.1
DO ialpha=2,nAlpha
  ALPHAS(ialpha) = 0.1 + ALPHAS(ialpha-1)
END DO



100 FORMAT('    ',E12.6,', ',E12.6,', ',E12.6,', ',E12.6,', ',E12.6,', ')

WRITE(10,*) "alphas = "
WRITE(10,100) ALPHAS
WRITE(10,*)
WRITE(10,*) "betas = "
WRITE(10,100) BETAS
WRITE(10,*)


WRITE(10,*) "S = "

DO ialpha=1,nAlpha
  ALPHA = ALPHAS(ialpha)
  PSQ = ALPHA*AM*temp

  DBW = EXP(-PSQ*GC(1))
  APS = PSQ*W1/AM
  BPS = PSQ

  DBWP = DBW/3.1416

  CALL SCINT(t,GC,GS,EPS,S1,temp, APS,BPS,DBWP,nTime,nBeta)

  SZCON = DBW*SQRT(AM/(12.566371*PSQ*W1*temp))
  DO i=1,nBeta
    S1(i)= S1(i)/EXP(BETAS(i)/2)
    S2(i) = SZCON*EXP(-AM*(EPS(i)**2+(PSQ*W1/AM)**2)/(4.*PSQ*temp*W1))
    S(i) = S1(i)+S2(i)
  ENDDO



  !if (ialpha.eq.1.or.ialpha.eq.2.or.ialpha.eq.20) then
  !  write(*,*)
  !  write(*,*) "Not Convolved", ALPHAS(ialpha)
  !  write(*,*) BETAS(1),S(1)
  !  write(*,*) BETAS(10),S(10)
  !  write(*,*) BETAS(26),S(26)
  !  write(*,*) BETAS(100),S(100)
  !endif


  
  DO i=1,nOsc
    RR = 0.5*X5(i)/temp
    U=EXP(RR)
    U=0.5*(U-1./U)
    ARG1(i)=W5*Q5(i)/(AM*X5(i)*U)
    ARG2(i)=W5*Q5(i)/(AM*X5(i)*TANH(RR))
    CALL BESSL(ARG1(i)*PSQ,BF,10, NMAX(i)) ! 10 first moments of modified 
    EX = EXP(-PSQ*ARG2(i))                 ! Bessel fct of first kind (I_n)
    DO j=1,10
      ANK(i,j)=BF(j)*EX
    ENDDO
  ENDDO

  CALL RCONV(nBeta,nOsc, NMAX,X5,ANK,temp,S1,BETAS)

  CALL ACON2(nBeta,NMAX,X5,ANK, temp, SZCON, EPS, AM, W1, PSQ, S2)
  
  DO i=1,nBeta
    S(i) = S1(i)+S2(i)
  ENDDO
  
  !write(*,*) ALPHA
  !DO i=1,nBeta
  !  WRITE(10,*) ALPHA, BETAS(i),S(i)
  !ENDDO
  WRITE(10,100) S



  !if (ialpha.eq.1.or.ialpha.eq.2.or.ialpha.eq.20) then
  !  write(*,*)
  !  write(*,*) "Convolved",ALPHAS(ialpha)
  !  write(*,*) BETAS(1),S(1)
  !  write(*,*) BETAS(10),S(10)
  !  write(*,*) BETAS(26),S(26)
  !  write(*,*) BETAS(100),S(100)
  !endif 
END DO

CLOSE(10)

END PROGRAM

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

SUBROUTINE ACON2(nBeta,NMAX,X5,ANK, temp, SZCON, EPS, AM, W1, PSQ, S2)

IMPLICIT NONE

integer :: i, j, k, l, m, ind1, ind2, nBeta, NMAX(2)
real(8) :: S2(nBeta),SK(1000), SINT, EIN, X5(2), SZCON, AM, PSQ, W1, temp, ANK(2,10)
real(8) :: EPS(nBeta)


  DO i=1,nBeta
    SK(i)=0.0
    l=NMAX(2)*2+1
    DO j=1,l
      m=NMAX(1)*2+1
      ind2 = j-NMAX(2)-1
      DO k=1,m
        ind1 = k-NMAX(1)-1
        EIN = ABS(EPS(i)-ind1*X5(1)-ind2*X5(2))
        SINT = SZCON*EXP(-AM*(EIN**2+(PSQ*W1/AM)**2)/(4.*PSQ*temp*W1))
        !if (SINT.GT.200) WRITE(*,*) SINT, ind1, ind2, ANK(2,ABS(ind2)+1), ANK(1,ABS(ind1)+1)
        SK(i) = SK(i) + ANK(2,ABS(ind2)+1)*ANK(1,ABS(ind1)+1)*SINT
      ENDDO
    ENDDO
    S2(i)=SK(i)
  ENDDO
END SUBROUTINE

SUBROUTINE RCONV(nBeta,nOsc, NMAX,X5,ANK,temp,S1,BETAS)

IMPLICIT NONE

INTEGER :: i, j, k, nBeta, nOsc, ind, l, NMAX(nOsc)

REAL(8) :: S1(nBeta), X5(nOsc), temp, BETAS(nBeta),  ANK(nOsc,20)
REAL(8) :: SLOG(1000), SK(1000), SINT, BETASIN

DO k=1,nOsc
  DO i=1,nBeta
    IF(S1(i).LE.0.0) THEN
      SLOG(i)=-100
    ELSE
      SLOG(i)=DLOG(S1(i))
    ENDIF
  ENDDO
  DO i=1,nBeta
    SK(i)=0.0
    l=NMAX(k)*2+1
    DO j=1,l
      ind = j-NMAX(k)-1
      BETASIN = ABS(BETAS(i)-FLOAT(ind)*X5(k)/temp)
      CALL STERP(BETASIN,BETAS,nBeta,SINT,SLOG)
      SK(i)=SK(i)+ANK(k,ABS(ind)+1)*SINT
    ENDDO
    S1(i)=SK(i)
  ENDDO
ENDDO

END SUBROUTINE


SUBROUTINE SCINT(t,GC,GS,EPS,S, temp, A, B ,F,NT,nBeta)

IMPLICIT NONE

INTEGER :: NT, nBeta
REAL(8) :: t(NT), GC(NT), GS(NT), S(nBeta),EPS(nBeta)
REAL(8) :: A, B, F, temp

INTEGER :: i, j

REAL(8) :: SM, SINSM, COSSM, V0, U, SINS, COSS, V, ST, CT, SINT, COST, AL
REAL(8) :: EX1, EX2
!REAL(8), ALLOCATABLE :: Q(:), R(:)
REAL(8) :: Q1,Q2,R1,R2

!ALLOCATE(Q(NT), R(NT))

EX1 = EXP(B*GC(1))
EX2 = EXP(-A*temp*t(1)**2)
!Q(1) = (COS(B*GS(1))*EX1-1.)*EX2
!R(1) = SIN(B*GS(1))*EX1*EX2
Q1 = (COS(B*GS(1))*EX1-1.)*EX2
R1 = SIN(B*GS(1))*EX1*EX2

DO i=1,nBeta
  AL = A-EPS(i)
  IF (AL.EQ.0) THEN
    WRITE(*,*) 'AL IS EQUAL to 0'
    !CALL INTG(t,Q,S(i),NT)
    !S(i)=S(i)*F
  ELSE
    S(i) = 0.
    SM = t(1)*AL
    SINSM = SIN(SM)
    COSSM = COS(SM)
    V0 = 0.
    DO j=2,NT
      U=t(j)*AL
      SINS = SIN(U)
      COSS = COS(U)
      V=U-SM
      IF (ABS(V/V0-1.).LE.5e-7) GO TO 40
      IF (ABS(V).LE.0.005) THEN
        ST = (V**2)/6.-(V**2)**2/120.
        CT = V*0.5-V**3/24.
      ELSE
        SINT = SINS*COSSM-COSS*SINSM
        COST = COSS*COSSM+SINS*SINSM
        ST = 1.-SINT/V
        CT = (1.-COST)/V
      ENDIF
  40  EX1 = EXP(B*GC(j))
      EX2 = EXP(-A*temp*t(j)**2)
!      Q(j) = (COS(B*GS(j))*EX1-1.)*EX2
!      R(j) = SIN(B*GS(j))*EX1*EX2
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

!DEALLOCATE(Q,R)

END SUBROUTINE



SUBROUTINE GTG(wgt,temp,AM,X,Q,t,PC,PS,nRho,NT,TBAR)

IMPLICIT NONE

INTEGER :: nRho, NT
REAL(8) :: X(nRho), Q(nRho), t(NT)
REAL(8) :: wgt, temp, AM, TBAR

INTEGER :: i
REAL(8) :: U, norm, F, H, A, C, CS, S
REAL(8) :: PC(NT), PS(NT)

U = Q(1)*X(1)/3.
CALL INTG(X,Q,A,nRho)
norm = wgt / AM / (U+A)
CALL FTRANS(temp,X,Q,t,PC,PS,nRho,NT)

F=X(1)*0.5/temp
CALL COTH(H,F)
DO i=2,NT
  U=X(1)*t(i)
  IF (U.LE.0.005) THEN
    C=0.5*U - U**3/24.
    S=U-U**3/6.
    CS = U/3. - U**3/30.
  ELSE
    C = COS(U)
    S = SIN(U)
    CS = S/U**2 - C/U
    C = (1.-C)/U
  ENDIF
  PC(i)=PC(i)+Q(1)/U*(H*(S-C)+C/F)
  PS(i)=PS(i)+Q(1)*CS
ENDDO
PC(1)= 0.5*Q(1)*(1./F+H)
PS(1)= 0.
TBAR = Q(1)*temp*X(1)/3.

DO i=1,nRho
  F=X(i)*0.5/temp
  CALL COTH(H,F)
  Q(i) = Q(i)*H/X(i)
ENDDO

CALL INTG(X,Q,A,nRho)
PC(1)=PC(1)+A

DO i=1,nRho
  Q(i) = Q(i)*X(i)**2*0.5
ENDDO
CALL INTG(X,Q,A,nRho)
TBAR= TBAR + A

DO i=1,NT
  PC(i) = PC(i)*norm
  PS(i) = PS(i)*norm
ENDDO

TBAR = TBAR*norm*AM/wgt

END SUBROUTINE

SUBROUTINE COTH(H,F)

IMPLICIT NONE

REAL(8) :: F, tmp, H

tmp = EXP(F)
H = (tmp+1./tmp)/(tmp-1./tmp)

END SUBROUTINE


SUBROUTINE INTG(X,Q,A,N)

! Trapeze integral

IMPLICIT NONE

INTEGER :: i, N

REAL(8) :: A
REAL(8) :: X(N),Q(N)

A=0.
DO i=2,N
  A = A + (Q(i)+Q(i-1))*(X(i)-X(i-1))*0.5
ENDDO

END SUBROUTINE

SUBROUTINE FTRANS(temp,X,Q,t,PC,PS,nRho,NT)

IMPLICIT NONE

INTEGER :: i,j, nRho, NT
REAL(8) :: PC(NT), PS(NT), t(NT)
REAL(8) :: X(nRho),Q(nRho),temp

REAL(8) :: S,SM, SINSM, COSSM, SINS, COSS, SINT, COST, ST, CT, Z, ZM, H, HM, U

DO i=2,NT
  PC(i)=0.
  PS(i)=0.
  SM=X(1)*t(i)
  SINSM = SIN(SM)
  COSSM = COS(SM)
  ZM = X(1)*0.5/temp
  DO j=2,nRho
    S = X(j)*t(i)
    SINS = SIN(S)
    COSS = COS(S)
    Z = X(j)*0.5/temp
    IF (ABS(S/SM-1.0).LE.5e-7) GO TO 40 !EXIT
    U = S-SM
    IF (U.LE.0.005) THEN
      ST = U**2/6. - U**4/120.
      CT = 0.5*U - U**3/24.
    ELSE
      SINT = SINS*COSSM - COSS*SINSM
      COST = COSS*COSSM + SINS*SINSM
      ST = 1.-SINT/U
      CT = (1.-COST)/U
    ENDIF
40  CALL COTH(H,Z)
    CALL COTH(HM,ZM)
    PC(i) = PC(i) + Q(j)/X(j)*H*(ST*SINS+CT*COSS) - Q(j-1)/X(j-1)*HM*(ST*SINSM-CT*COSSM)
    PS(i) = PS(i) + Q(j)/X(j)*(CT*SINS-ST*COSS) + Q(j-1)/X(j-1)*(CT*SINSM+ST*COSSM)
    SM = S
    SINSM=SINS
    COSSM = COSS
    ZM=Z
  ENDDO
  PC(i)=PC(i)/t(i)
  PS(i)=PS(i)/t(i)
ENDDO

END SUBROUTINE

SUBROUTINE STERP(B,BETAS,NB,SINT, SLOG)

IMPLICIT NONE

INTEGER :: IC, NB
REAL(8):: B, BETAS(NB), SLOG(150), SINT, SL

IF(B.LE.BETAS(1)) THEN
  SINT = EXP(SLOG(1))
ELSE IF (B.GE.BETAS(NB)) THEN
  SINT = 0.
ELSE
  IC = 1
10  IF (B.GE.BETAS(IC)) GOTO 20
  IC = IC -1
  GOTO 10
20  IF (B.LT.BETAS(IC+1)) GOTO 30
  IC = IC +1
  GOTO 10
30 SL = SLOG(IC)+((B-BETAS(IC))/(BETAS(IC+1)-BETAS(IC)))*(SLOG(IC+1)-SLOG(IC))
  SINT=EXP(SL)
ENDIF

END SUBROUTINE

SUBROUTINE BESSL(X,B,NX,NMAX)

IMPLICIT NONE

REAL(8) :: X, B(10), BF(60), TA, FNFACT, FACT, X2N, X2
INTEGER :: NX, IND, IOTA, IORD, i, j, NMAX

NMAX = 0.0

IF (X.LT.0.05) THEN
  FNFACT=1.0
  X2N=1.0
  X2=X*0.5
  DO i=1,NX
    B(i)=X2N/FNFACT
    X2N=X2N*X2
    FNFACT=FNFACT*FLOAT(i)
    IF(B(i).LT.1.e-20) THEN
      B(i) = 0.0
      IF (NMAX.EQ.0) NMAX = i-1
    ENDIF
  ENDDO
ELSE
  DO i=1,NX
    B(i)=0.0
    IF((X-1.0).LT.0) THEN
      IORD=-37.0/(.43429*DLOG(.1*X))
      IF((IORD-5).LT.0) B(1)=1
    ELSE
      IORD=30
    ENDIF
    BF(IORD-1)=1.e-37
    BF(IORD)=0.0
    TA=1.e-37
    IOTA = IORD-1
    DO j=2,IOTA
      IND = IORD-j
      BF(IND)=FLOAT(IND)*2.0*BF(IND+1)/X+BF(IND+2)
      TA=TA+BF(IND)
    ENDDO
    TA=2.0*TA-BF(1)
    FACT=EXP(X)/TA
    DO j=1,NX
      B(j)=FACT*BF(j)
      IF(B(j).LT.1.e-20) THEN
        B(j) = 0.0
        IF (NMAX.EQ.0) NMAX = j-1
      ENDIF
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE
