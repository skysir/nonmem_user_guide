$PROBLEM WARFARIN
$INPUT ID TIME AMT RATE EVID MDV DV TSTRAT TMIN TMAX
$DATA warfarin_pfim.csv ignore=C

$SUBROUTINES ADVAN2 TRANS1

$PK
KA=THETA(1)*EXP(ETA(1))
K=THETA(2)*EXP(ETA(2))
V=THETA(3)*EXP(ETA(3))
S2=V
F1=1.0

$ERROR
IPRED=A(2)/V
;W=THETA(4)+THETA(5)*IPRED
;Y=IPRED + W*EPS(1)
Y=IPRED + EPS(1) + IPRED*EPS(2)

$THETA
2.0 ;[KA]
0.25  ;[K]
15.0  ;[V]
;0.5  ; [siga]
;0.15 ; [SIGP]

$OMEGA (1.0) (0.25) (0.10)
;$SIGMA (1.0 FIXED)
$SIGMA 0.25 0.0225

$DESIGN GROUPSIZE=200 FIMTYPE=1 VARCROSS=1 APPROX=FO
