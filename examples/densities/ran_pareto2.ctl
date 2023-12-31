$PROB TESTING RANDOM SAMPLERS
$ABBR FUNCTION pareto2(VQI,10)
$ABBR VECTOR VV(10)
$INPUT  AMT TVAL DV
$DATA  rsampler.csv

$PRED
YMIN=theta(1)
lambda=theta(2)
ALPHA=THETA(3)
IF(ICALL==4) THEN
VV(1)=1.0
VV(2)=YMIN
VV(3)=LAMBDA
VV(4)=ALPHA
" CALL pareto2_RNG(2,VV)
Y=VV(1)
DV=Y
ELSE
VQI(1)=DV
VQI(2)=YMIN
VQI(3)=LAMBDA
VQI(4)=ALPHA
WW=pareto2(VQI)
Y=2.0*WW
ENDIF

$THETA (3.0 FIXED) (0.0,0.6) (0.0,1.3)

$SIMULATION (567811 NORMAL) (2933012 UNIFORM) (445678 NORMAL) SUBPROBLEMS=1 
$EST METHOD=0 MAXEVAL=9999 PRINT=1 -2LL NOTHETABOUNDTEST NOABORT
$COVR
$TABLE TVAL DV NOAPPEND NOPRINT FILE=ran_pareto2.tab
