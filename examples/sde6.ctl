; Control problem for sde system.  No SDE matters done.  From Chris Tornoe
$PROBLEM PK ODE HANDS ON ONE

$INPUT ID TIME DV AMT CMT FLAG

$DATA   sde6.csv
        IGNORE=@

$SUBROUTINE ADVAN6 TOL 10 DP

$MODEL 
       COMP = (CENTRAL);

$PK
  IF(NEWIND.NE.2) OT = 0
   
  TVCL  = THETA(1)
  CL    = TVCL*EXP(ETA(1))
  
  TVVD  = THETA(2)
  VD    = TVVD*EXP(ETA(2))


$DES
 DADT(1) = - CL/VD*A(1) ;+SGW1
 
$ERROR 
  
     IPRED = A(1)/VD
     IRES  = DV - IPRED
     W     = THETA(3)
     IWRES = IRES/W
     Y     = IPRED+W*EPS(1)

$THETA (0,10)               ;1 CL
$THETA (0,32)               ;2 VD
$THETA (0, 2)               ;4 SIGMA

$OMEGA 0.1                  ;1 CL
$OMEGA 0.01                 ;2 VD

$SIGMA 1 FIX                ; PK

$EST MAXEVAL=9999 METHOD=1 LAPLACE NUMERICAL SLOW INTER NOABORT SIGDIGITS=3 PRINT=1 MSFO=sde6.msf
$COV MATRIX=R

$TABLE ID TIME FLAG AMT CMT IPRED IRES IWRES
       ONEHEADER NOPRINT FILE=sde6.fit
