$PROB RUN# 
$INPUT C ID TIME DV AMT RATE EVID MDV CMT ROWNUM SID
$DATA superid3_sim.csv IGNORE=C

$SUBROUTINES ADVAN2 TRANS2

$PK
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
KA=DEXP(MU_1+ETA(1)+ETA(4))
CL=DEXP(MU_2+ETA(2)+ETA(5))
V=DEXP(MU_3+ETA(3)+ETA(6))
S2=V

$ERROR
IPRE=F
Y = IPRE + IPRE*EPS(1)

; Initial values of THETA
$THETA 0.18 -5.3 -3.0
;INITIAL values of OMEGA
$OMEGA BLOCK(3)
0.01
0.001 0.01
0.001 0.001 0.01

$OMEGA BLOCK(3)
0.03
0.001 0.03
0.001 0.001 0.03

;Initial value of SIGMA
$SIGMA 
0.003     ;[P]

$LEVEL
SID=(4[1],5[2],6[3])

$SIMULATION (567811 NORMAL) (2933012 UNIFORM) ONLYSIMULATION SUBPROBLEMS=1
$TABLE C ID TIME DV AMT RATE EVID MDV CMT ROWNUM SID KA CL V ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
NOAPPEND ONEHEADER FILE=superid3.csv  NOPRINT
      