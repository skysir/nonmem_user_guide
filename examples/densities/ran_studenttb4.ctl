$PROB TESTING RANDOM SAMPLERS
$ABBR FUNCTION STUDENTTB2(VQI,10)
$ABBR VECTOR VV(10)
$INPUT  AMT TVAL DV DV2
$DATA  rsampler.csv

$PRED
NU=THETA(1)
QM1=THETA(2)
QM2=THETA(3)
SC11=THETA(4)
SC12=THETA(5)
SC22=THETA(6)
IPRED=QM1
F=QM1
IPRED2=QM2
IF(ICALL==4) THEN
VV(1)=1.0
VV(2)=1.0
VV(3)=NU
VV(4)=QM1
VV(5)=QM2
VV(6)=SC11
VV(7)=SC12
VV(8)=SC22
; DESCRIPTION : Given Normal Random generator K, and parameters X(3)..., 
;               return two random univariate t-samples X(1),X(2) followed by superimposed correlation C12=X(7)
" CALL STUDENTTB4_RNG(3,VV)
Y=VV(1)
DV=Y
DV2=VV(2)
ELSE
VQI(1)=DV
VQI(2)=DV2
VQI(3)=NU
VQI(4)=QM1
VQI(5)=QM2
VQI(6)=SC11
VQI(7)=SC12
VQI(8)=SC22
WW=STUDENTTB2(VQI)
Y=2.0*WW
ENDIF

$THETA 4.0 30.0 35.0 (0.0,10.0) (0.8) (0.0,12.0)

$SIMULATION (567811 NORMAL) (2933012 UNIFORM) (445678 NORMAL) SUBPROBLEMS=1
$EST METHOD=0 MAXEVAL=9999 PRINT=10 -2LL NOHABORT
$COVR
$TABLE TVAL DV DV2 IPRED IPRED2 NOAPPEND NOPRINT FILE=ran_studenttb4.tab
