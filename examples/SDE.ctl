$PROBLEM ODE MODEL

$INPUT ID TIME DV AMT CMT FLAG SDE MDV

$DATA   sde.csv
        IGNORE=@

;$MSFI runODE.msf NPOPETAS=3

$SUBROUTINE ADVAN13 TOL 6 OTHER=sde.f90

$MODEL COMP = (ABS);
       COMP = (CENTRAL);
       COMP = (DFDX1)
       COMP = (DFDX2)
       COMP = (DPDT11)
       COMP = (DPDT21)
       COMP = (DPDT22)
$PK
   
  TVCL  = THETA(1)
  CL    = TVCL*EXP(ETA(1))
  
  TVVD  = THETA(2)
  VD    = TVVD*EXP(ETA(2))
  
  THL  = THETA(3)
  HLA  = THL*EXP(ETA(3))
  
  KA = LOG(2)/HLA
  SGW3=THETA(5)

$DES 
" FIRST
" REAL*8 SGW(2)
 DADT(1) = - KA*A(1)
 DADT(2) = KA*A(1) - CL/VD*A(2)
; NEXT DERIVATIVES ARE ACUALLY PREDICTIVE VALUES FOR COMPARTMENTS 1 AND 2, RESPECTIVELY
;  Derivatives of these with respect to A() will be calcualted symbolically by DES routine created by NMTRAN
 DADT(3) = A(1)
 DADT(4) = LOG((DABS(A(2))+1.0D-300)/VD)
; DUMMY PLACEMENT FOR DERIVATIVES OF THE STOCHASTIC ERROR SYSTEM.  THESE ARE FILLED OUT BY SDE_DER
" SGW(1)=0.0D+00
" SGW(2)=SGW3
;  the DA() array THEN contains all derivatives of DADT (=DXDT) with respect to A(=X).
; Number of compartments=2, number of base model derivative equations=2
"LAST
"      CALL SDE_DER(DADT,A,DA,IR,SGW,2.0d+00,2d+00)
 
$ERROR
     IF(ICALL.EQ.4) THEN
       IF(DV.NE.0) Y = LOG(DV)
       RETURN
     ENDIF
     IPRED = LOG((DABS(A(2))+1.0D-300)/VD)
     W     = THETA(4)
     IRES  = DV - IPRED
     IWRES = IRES/W
     WS=0.0
     IF(CMT==2 .OR. CMT==0) THEN
; CENTRAL COMPARTMENT, PLASMA LEVELS
; EPS(2) = USER MODEL ERROR CONTRIBUTION
; EPS(4) = STOCHASTIC ERROR CONTRIBUTION.  THE WS IS JUST A PLACEHOLDER COEFFICIENT.  SDE_CADD WILL REPLACE THIS
; WITH THE CORRECT VALUE
     Y     = IPRED+W*EPS(2) + WS*EPS(4)
     ENDIF
     IF(CMT==1 ) THEN
; ABSORPTION COMPARTMENT.  IN THIS PROBLEM, THERE ARE NO OBSERVED VALUES FOR ABSORPTION COMPARTMENT.
; BUT PUT THESE IN AS PLACE HOLDERS ANYWAY.
; EPS(1) = USER MODEL ERROR CONTRIBUTION
; EPS(3) = STOCHASTIC ERROR CONTRIBUTION
     Y     = IPRED+W*EPS(1) + WS*EPS(3)
     ENDIF
; SDE_CADD WILL EVALUATE THE TRUE COEFFICIENTS (WS) TO THE STOCHASTIC COMPONENTS.
; Number of compartments=2, number of base model derivative equations=2
"LAST
"       CALL SDE_CADD(A,HH,TIME,DV,CMT,2.0D+00,2.0D+00,SDE)

$THETA (0,22.0    )     ; 1  CL
$THETA (0,62.0    )     ; 2  VD
$THETA (0,35.0    )     ; 3  HLA
$THETA (0, 0.2     )     ; 4  PK PROP MEAS ERROR
$THETA (0, 1.0)          ; 5 SGW3

$OMEGA 0.1               ;1 CL
$OMEGA 0.1                ;2 VD
$OMEGA 0.1                ;3 HLA

$SIGMA 1 FIX               ; PK-ABSORPTION
$SIGMA 1 FIX               ; PK-CENTRAL
$SIGMA 1 FIX               ; PK-ABSORPTION, STOCHASTIC.
$SIGMA 1 FIX               ; PK-CENTRAL, STOCHASTIC.

$SIM(1)
$EST MAXEVAL=9999 METHOD=1 INTERACTION SIGDIGITS=3 PRINT=1 NOABORT MSFO=runODE.msf
;$COV

;$TABLE ID TIME FLAG AMT CMT IPRED IRES
;       ONEHEADER NOPRINT FILE=sdtabODE

;$TABLE ID CL VD HLA ETA1 ETA2 ETA3
;       ONEHEADER NOPRINT FILE=patabODE