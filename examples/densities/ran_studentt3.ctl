$PROB TESTING RANDOM SAMPLERS
$ABBR FUNCTION STUDENTT(VQI,5)
$ABBR VECTOR VV(5)
$INPUT  AMT TVAL DV
$DATA  rsampler.csv

$PRED
NU=THETA(1)
QM=THETA(2)
SIGV=THETA(3)
IPRED=QM
F=QM
IF(ICALL==4) THEN
VV(1)=1.0
VV(2)=NU
VV(3)=QM
VV(4)=SIGV
; Create t with a normal deviate divided by square root of chi-square deviate
" CALL STUDENTT3_RNG(3,2,VV)
Y=VV(1)
DV=Y
ELSE
VQI(1)=DV
VQI(2)=NU
VQI(3)=QM
VQI(4)=SIGV
VQI(5)=1.0
WW=STUDENTT(VQI)
Y=2.0*WW
ENDIF

$THETA 4.5 30.0  (0.0,10.0)
;$OMEGA (0.0 fixed)

$SIMULATION (567811 NORMAL) (2933012 UNIFORM) (445678 NORMAL) SUBPROBLEMS=1
$EST METHOD=0 MAXEVAL=9999 PRINT=10 -2LL NOHABORT
$COVR
$TABLE TVAL DV IPRED NOAPPEND NOPRINT FILE=ran_studentt3.tab
