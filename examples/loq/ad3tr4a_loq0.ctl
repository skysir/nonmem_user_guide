$PROB  AD3TR4
$INPUT  C SET ID JID TIME CONC AMT RATE EVID MDV CMT DV LOQ TYPE

$DATA  ad3tr4a.csv IGNORE = C
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
LAQ=3.0
SD = THETA(5)
IPRED = LOG(F)
DUM = (LOQ - IPRED) / SD
CUMD = PHI(DUM)+1.0E-10
DUMA = (LAQ - IPRED) / SD
CUMDA = PHI(DUMA)-1.0E-10
IF(TYPE.EQ.2) DV_LOQ=LOQ
IF(TYPE.EQ.3) DV_LAQ=LAQ
IF (TYPE .EQ. 1.OR.NPDE_MODE==1) THEN
      F_FLAG = 0
      Y = IPRED + SD * ERR(1)
ENDIF
IF (TYPE .EQ. 2.AND.NPDE_MODE==0) THEN
      F_FLAG = 1
      Y = CUMD
      MDVRES=1
ENDIF
IF (TYPE .EQ. 3.AND.NPDE_MODE==0) THEN
      F_FLAG = 1
      Y = (1.0-CUMDA)
      MDVRES=1      
ENDIF
;IF(TYPE==3.0.AND.NPDE_MODE==1) WRITE(*,*) CUMDA

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

;Initial value of SIGMA
$SIGMA 
1.0 FIXED   ;[P]

$EST METHOD=COND INTERACTION LAPLACE MAXEVAL=9999 NSIG=3 SIGL=10 PRINT=5 NOHABORT MCETA=10
;$COV MATRIX=R PRINT=E UNCONDITIONAL
$TABLE ID TIME DV IPRED PRED RES EPRED EIPRED EIRES EIWRES CPRED CWRES ECWRES EWRES TYPE NPDE NPD CIWRES NOAPPEND ONEHEADER 
 ESAMPLE=1000 NPDTYPE=0 FILE=ad3tr4a_loq0.TAB NOPRINT SEED=16993234
