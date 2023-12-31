$PROBLEM PDLIDR
$INPUT ID AMT TIME DV EVID  CMT
$DATA PDLIDR.csv IGNORE=C
$SUBROUTINES ADVAN16 TOL=4
$MODEL NCOMPARTMENTS=3 

$PK
MU_1=LOG(THETA(3))
MU_2=LOG(THETA(4))
MU_3=LOG(THETA(5))
MU_4=LOG(THETA(6))
MU_5=LOG(THETA(7))
MXSTEP=2000000000
KEL=THETA(1)
V=THETA(2)
K0=EXP(MU_1+ETA(1))
K1=EXP(MU_2+ETA(2))
SMAX=EXP(MU_3+ETA(3))
SC50=EXP(MU_4+ETA(4))
;  TAUy
TAU1=EXP(MU_5+ETA(5))
; Initial conditions
A_0(1)=0
A_0(2)=K0/K1
A_0(3)=K0*TAU1

$DES
; AD_x_y is the State value of A(x) delayed for time TAUy.  
; AP_x_y is the State value of A(x) in the past, for time delay TAUy.  
AP_2_1=K0/K1
;BASE EQUATIONS 
CC=A(1)/V
DADT(1)=-KEL*A(1)
DADT(2)=K0*(1+SMAX*CC/(SC50+CC))-K1*A(2)
DADT(3)=K1*A(2)-K1*AD_2_1

$ERROR

Y1=LOG(A(3))  
IF(CMT==3) IPRED=Y1
IF(CMT==3) Y=IPRED+EPS(1)

$THETA
0.25 FIX    ; 1: KEL
1    FIX    ; 2: V
(0,0.5,5)   ; 3: K0
(0,0.05,0.5); 4: K1
(0,50,500)  ; 5: SMAX
(0,1,10)    ; 6: SC50
(5,20,200)  ; 7: TR

$OMEGA BLOCK(5) VALUES(0.2,0.001)

$SIGMA
0.01
$EST METHOD=SAEM INTERACTION ISAMPLE=2 NBURN=200 NITER=200 PRINT=10 NOHABORT CTYPE=3 SIGL=3 RANMETHOD=3S2P
$EST METHOD=IMP INTERACTION MAPITER=0 NITER=200 ISAMPLE=300 CTYPE=3 PRINT=1
$COV UNCONDITIONAL MATRIX=R
$TABLE ID TIME IPRED EVID CMT NOPRINT ONEHEADER
FILE=PDLIDR.tab
