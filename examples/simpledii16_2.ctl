;$SIZES MAXNRDS=1 PAST_SIZE=10000
$PROBLEM RA
$ABBR DERIV2=NO DERIV2=NOCOMMON DES=FULL ; DERIV1=NO
$INPUT ID TIME  AMT  RATE II ADDL CMT  EVID MDV  DV
$DATA simpleDii16_2.csv IGNORE=C
$SUBROUTINES ADVAN16 TOL=10 ATOL=10
$MODEL NCOMPARTMENTS=4

$PK
K10=THETA(1)+ETA(1)
K12=THETA(2)+ETA(2)
K21=THETA(3)+ETA(3)
V1=1.0
S1=V1
KEO=THETA(4)
K40=THETA(5)
TAU1=THETA(6)
" AQ=0.0
	
$DES
MXSTEP=2000000000
PASTZERO=-240.0
 AP_1_1=0.0

 DADT(1)=-K10*A(1)-K12*A(1)+K21*A(2) 
 DADT(2)=K12*A(1)-K21*A(2)
 DADT(3)=KEO*A(1)-KEO*AD_1_1
 DADT(4)=KEO*AD_1_1-K40*A(4)

 AQ=AD_1_1

$ERROR

IPRED=A(1)/V1
Y=IPRED*(1.0+EPS(1))
A1=A(1)
A2=A(2)
A3=A(3)
A4=A(4)
A5=0.0
" A5=AQ

$THETA
0.17 ; [K10]
0.37 ; [K12]
0.25 ; [K21]
0.2   ; KEO
0.05  ; K40
4.0 ; [TAU1]

$OMEGA (0.0 FIXED) (0.0 FIXED) (0.0 FIXED)
$SIGMA (0.0 FIXED)

$SIMULATION (567811 NORMAL) ONLYSIMULATION SUBPROBLEMS=1

$TABLE ID TIME A1 A2 A3 A4 A5 EVID NOAPPEND NOPRINT ONEHEADER
FILE=simpledii16_2.tab
