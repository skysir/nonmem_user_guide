$PROBLEM WARFARIN
$INPUT ID TIME AMT RATE EVID MDV DV TSTRAT TMIN TMAX REPX
$DATA warfarin.csv ignore=C ; REPL=5

$SUBROUTINES ADVAN2 TRANS2

$PK
MU_1=LOG(THETA(1))
MU_2=LOG(THETA(2))
MU_3=LOG(THETA(3))
CL=THETA(1)*EXP(ETA(1))
V=THETA(2)*EXP(ETA(2))
KA=THETA(3)*EXP(ETA(3))
S2=V
F1=1.0

$ERROR
IPRED=A(2)/V
Y=IPRED*(1.0+EPS(1))

$THETA
0.15 ;[CL]
8.0  ;[V]
1.0  ;[KA]

$OMEGA (0.07 ) (0.02) (0.6 )
$SIGMA (0.01 )

;$OMEGA (1.61) (0.46) (13.8)
;$SIGMA (5.714)

$OPTDESIGN MODE=0 FIMDIAG=1 OFVTYPE=0 GROUPSIZE=1.0 NELDER DESEL=TIME DESELSTRAT=TSTRAT DESELMIN=TMIN DESELMAX=TMAX SEED=4455322
$EST METHOD=0 MAXEVAL=0 SIGL=12 nohabort PRINT=10 ; FORMAT=S1PE23.16
$COV MATRIX=R UNCONDITIONAL SIGL=12 CHOLROFF=0
$TABLE ID TIME EVID MDV DV IPRED NOPRINT NOAPPEND VARCALC=3 FILE=warfarin.tab   ; FORMAT=S1PE23.16
