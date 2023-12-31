$PROB  AD3TR4
$INPUT  C SET ID JID TIME CONC AMT RATE EVID MDV CMT DV LOQ TYPE

$DATA  ad3tr4.csv IGNORE = C ; IGNORE=(TYPE=2)
$SUBROUTINES ADVAN3 TRANS4

$PK
MU_1 = THETA(1)
MU_2 = THETA(2)
MU_3 = THETA(3)
MU_4 = THETA(4)
CL = EXP(MU_1 + ETA(1))
V1 = EXP(MU_2 + ETA(2))
Q = EXP(MU_3 + ETA(3))
V2 = EXP(MU_4 + ETA(4))
S1=V1

$ERROR
SD = THETA(5)
IPRED = LOG(F)
DUM2 = (DV - IPRED) / SD
DUM = (LOQ - IPRED) / SD
CUMD = PHI(DUM)+1.0E-30
CUMD2 = PHI(DUM2)+1.0E-30
IF(TYPE.EQ.1) THEN
Y=2.0*LOG(SD)+DUM2*DUM2
CDF_L=CUMD2
ENDIF
IF(TYPE.EQ.2) THEN
Y = -2.0*LOG(CUMD)
CDF_L=CUMD
DV_LOQ=LOQ
ENDIF

$THETA 
         1.67E+00  1.60E+00  7.76E-01  2.36E+00  2.75E-01

;INITIAL values of OMEGA
$OMEGA BLOCK(4)
0.15   ;[P]
0.01  ;[F]
0.15   ;[P]
0.01  ;[F]
0.01  ;[F]
0.15   ;[P]
0.01  ;[F]
0.01  ;[F]
0.01  ;[F]
0.15   ;[P]

$EST METHOD=COND LAPLACE -2LL MAXEVAL=9999 NSIG=3 SIGL=10 PRINT=5 NOHABORT MCETA=10
;$COV MATRIX=R PRINT=E UNCONDITIONAL
$TABLE ID TIME DV IPRED NPD NOAPPEND ONEHEADER ESAMPLE=1000
 FILE=ad3tr4_loq6.TAB NOPRINT
