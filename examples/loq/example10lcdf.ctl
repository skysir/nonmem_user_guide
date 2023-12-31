$PROB  F_FLAG04est2a.ctl
$INPUT C ID DOSE=AMT TIME DV WT TYPE
$DATA example10lcdf.csv IGNORE=@

$SUBROUTINES  ADVAN2 TRANS2


$PK
   CALLFL=1
   MU_1=THETA(1)
   KA=DEXP(MU_1+ETA(1))
   MU_2=THETA(2)
   V=DEXP(MU_2+ETA(2))
   MU_3=THETA(3)
   CL=DEXP(MU_3+ETA(3))
   SC=V/1000

$THETA  1.6 2.3 0.7 0.1 0.1

$OMEGA BLOCK (3)
0.5
0.01 0.5
0.01 0.01 0.5


; Because THETA(4) and THETA(5) have no inter-subject variability 
; associated with them, the algorithm must use a more computationally 
; expensive gradient evaluation for these two parameters

$SIGMA 0.1


$PRIOR NWPRI
; Priors to Omegas
$OMEGAP BLOCK (3)
0.09 FIX
0.0 0.09
0.0 0.0 0.09
$OMEGAPD (3 FIX)


$ERROR
EXCL2=1.0-TYPE
EXCL=TYPE
EXCL3=0.0
IF(EVID/=0) EXCL=1.0
IF(EVID/=0) EXCL2=1.0
IF(EVID/=0) EXCL3=1.0
    EXPP=THETA(4)+F*THETA(5)
IPRED=F
; Use protected exponent PEXP, to avoid numerical overflow
A=PEXP(EXPP)
B=1.0+A
IF (TYPE.EQ.0.OR.NPDE_MODE==1) THEN
; PK Data
    F_FLAG=0
    Y=F+F*ERR(1) ; a prediction
 ELSE
; Categorical data
    F_FLAG=1
    Y=DV*A/B+(1.0-DV)/B      ; a likelihood
    MDVRES=1
 ENDIF
IF(TYPE==1)  THEN
CDF_L=(1.0-DV)*1.0/B + DV
CDF_LA=DV*1.0/B
DV_LOQ=DV
DV_LAQ=DV-1.0
ENDIF

$EST METHOD=ITS INTER LAP NITER=30 PRINT=5 SIGL=6 NSIG=2 NOHABORT 
     NOPRIOR=1 CTYPE=3 CITER=10 CALPHA=0.05
;$EST METHOD=COND LAP INTER MAXEVAL=9999 PRINT=1 NOPRIOR=1
$COV UNCONDITIONAL PRINT=E MATRIX=R SIGL=10
$TABLE ID TIME DV IPRED PRED RES EPRED EIPRED EIRES EIWRES CPRED CWRES ECWRES EWRES TYPE NPDE NPD CIWRES NOAPPEND ONEHEADER 
 ESAMPLE=1000 NPDTYPE=0 EXCLUDE_BY EXCL3 FILE=example10lcdf.TAB NOPRINT SEED=16993234
$TABLE ID TIME DV IPRED PRED RES EPRED EIPRED EIRES EIWRES CPRED CWRES ECWRES EWRES TYPE NPDE NPD CIWRES NOAPPEND ONEHEADER 
 ESAMPLE=1000 NPDTYPE=0 EXCLUDE_BY EXCL FILE=example10lcdf0.TAB NOPRINT SEED=16993234
$TABLE ID TIME DV IPRED PRED RES EPRED EIPRED EIRES EIWRES CPRED CWRES ECWRES EWRES TYPE NPDE NPD CIWRES NOAPPEND ONEHEADER 
 ESAMPLE=1000 NPDTYPE=0 EXCLUDE_BY EXCL2 FILE=example10lcdf1.TAB NOPRINT SEED=16993234
