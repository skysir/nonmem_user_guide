$PROB  F_FLAG04est2a.ctl
$INPUT C ID DOSE=AMT TIME DV WT TYPE
$DATA example10.csv IGNORE=@

$SUBROUTINES  ADVAN2 TRANS2


$PK
   CALLFL=1
   MU_1=DLOG(THETA(1))
   KA=DEXP(MU_1+ETA(1))
   MU_2=DLOG(THETA(2))
   V=DEXP(MU_2+ETA(2))
   MU_3=DLOG(THETA(3))
   CL=DEXP(MU_3+ETA(3))
   SC=V/1000

$THETA  5.0 10.0 2.0 0.1 0.1

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
    EXPP=THETA(4)+F*THETA(5)
IF (TYPE.EQ.0) THEN
; PK Data
    F_FLAG=0
    Y=F+F*ERR(1) ; a prediction
 ELSE
; Categorical data
    F_FLAG=1
; Use protected exponent PEXP, to avoid numerical overflow
    A=PEXP(EXPP)
    B=1+A
    Y=DV*A/B+(1-DV)/B      ; a likelihood
 ENDIF



$EST METHOD=ITS INTER LAP NITER=1000 PRINT=5 SIGL=6 NSIG=2 
     NOABORT NOPRIOR=1 CTYPE=3 CITER=10 CALPHA=0.05 
     FILE=example10.ext
; Because of categorical data, which can make conditional density highly 
; non-normal, select a t-distribution with 4 degrees of freedom for 
; importance sampling proposal density
$EST METHOD=IMP INTER LAP NITER=1000 PRINT=1 ISAMPLE=300 DF=4 
     IACCEPT=1.0
$EST METHOD=IMP EONLY=1 NITER=5 ISAMPLE=1000 PRINT=1 DF=4 
     IACCEPT=1.0 MAPITER=0 

$EST METHOD=SAEM EONLY=0 INTER LAP NBURN=2000 NITER=1000 PRINT=50 
     DF=0 IACCEPT=0.4
$EST METHOD=IMP EONLY=1 NITER=5 ISAMPLE=1000 PRINT=1 DF=4 
     IACCEPT=1.0 MAPITER=0 

$EST METHOD=BAYES NBURN=3000 NSAMPLE=3000 PRINT=100 
     FILE=example10.txt DF=0 IACCEPT=0.4 NOPRIOR=0

$EST METHOD=COND LAP INTER MAXEVAL=9999 PRINT=1 FILE=example10.ext
     NOPRIOR=1 NOHABORT

$COV UNCONDITIONAL PRINT=E MATRIX=R SIGL=10
$TABLE ID DOSE WT TIME TYPE DV A NOPRINT FILE=example10.tab
