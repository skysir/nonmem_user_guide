;DDE
$PROBLEM LOGISTIC
; turn off second derivative assessments, sometimes even 1st derivatives if only simulating
$ABBR DERIV2=NO DERIV2=NOCOMMON 
$INPUT ID AMT TIME PRDV DV EVID MDV 
$DATA LOGISTIC6.csv IGNORE=C
$SUBROUTINES ADVAN16 TOL=6 ATOL=6 
$MODEL NCOMPARTMENTS=1

$PK
CALLFL=-2
MXSTEP=2000000000
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
KG=EXP(MU_1+ETA(1))
Y0=EXP(MU_2+ETA(2))
YSS=EXP(MU_3+ETA(3))
TAU1=EXP(MU_4+ETA(4))
; Initial conditions
A_0(1)=Y0



TSTOP=500.0
; INITIALIZING EQUATIONS FOR DDE COMPARTMENTS
TAU_1=1*TAU1

$DES
; AD_1_1 is the State value of A(1) delayed for time TAU1.
; AP_1_1 is the State value of A(1) in the past, for time delay TAU1.

; DELAY SETUP FOR EQUATION SET 1
 AP_1_1=Y0
; DELAY EQUATIONS FOR EQUATION SET 0 (BASE EQUATIONS)
;BASE EQUATIONS
 DADT(1)=KG*(1.0-AD_1_1/YSS)*A(1)

$ERROR
A1=A(1)


Y1=1.0
IPRED=A(1)
Y=IPRED*(1.0+EPS(1))

;$THETA
;-1.609     ; KG
;-0.0001     ; Y0
;2.3026    ; YSS
;1.609     ; TAU1

;$OMEGA (0.01)x4

;$SIGMA
;0.003


$THETA
-1.2
-2.34841E-02 
2.8
0.9

$OMEGA BLOCK(4)
0.1
6.83328E-04  0.1
1.05324E-03 -1.82494E-03  0.1
7.11753E-04  8.46056E-04  1.47772E-03  0.1

$SIGMA
0.01

$EST METHOD=ITS INTERACTION NOABORT SIGL=4 SIGLO=6 MCETA=10 NSIG=2 PRINT=1 NITER=100 CTYPE=3 FAST
$EST METHOD=IMP INTERACTION MAXEVAL=9999 NOABORT SIGL=6 NSIG=2 PRINT=1 NITER=100 CTYPE=3 MAPITER=0
$COV MATRIX=R UNCONDITIONAL
