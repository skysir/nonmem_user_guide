$PROB BIVARIATE EXAMPLE
; THESE DECLARATIONS ALLOW ANY FUNCTION TO HAVE ALTERNATIVE DIMENSIONS FOR THEIR ARRAYS
; BUT, USER DEFINED DIMENSIONS ARE PASSED AS THE LAST ARGUMENT TO FUNC, SUCH AS:
; BV=BIVARIATE(VBI(1),FNC001_1(1,1),FNC001_2(1,1,1),5)
$ABBR FUNCTION BIVARIATE(VBI,5)

$INPUT SIM ID DOSE DV TIME
$DATA bivariate.csv IGNORE=C

$SUBROUTINES OTHER=bivariate.f90

$PRED
  B1=THETA(1)
  B2=THETA(2)
  B3=THETA(3)
  K =LOG(2)/EXP(THETA(4))
  ED50=EXP(THETA(5))
  U =(1-EXP(-K*TIME ))

  MU_1=B1+B3*DOSE/(DOSE+ED50)
  MU_2=B2
  MXB=MU_1+ETA(1)
  MXU=MU_2+ETA(2)
  MX =MXB + MXU*U  ;***Current model prediction***;

  PHIMX=PHI(MX)

   IF(NEWIND.NE.2) THEN
   TIMEP=0
   MXP=0
   DVP=0
   PHIMXP=0.5
   ENDIF

  RHOB=(2/(1+EXP(-THETA(8)))-1)
  IF(RHOB>0.0)  RHO=RHOB**(TIME-TIMEP)
  IF(RHOB==0.0) RHO=0.0
  IF(RHOB<0.0)  RHO=-(-RHOB)**(TIME-TIMEP)

  PC =(1-PHIMX) *(1-DV ) + PHIMX*DV
  IF(PC.LE.0.0) EXIT


  V=SQRT(1+OMEGA(1,1)+OMEGA(2,2)*U**2)
  POPP = (B1+B2*U +B3*DOSE/(DOSE+ED50))/V ;***Population mean prediction***;

  IF (TIME.EQ.1) THEN
    JP=PC
    PCP=1.0
  ELSE
  ;***Pass information to bivariate normal***;
   VBI(1)=RHO
   VBI(2)=MX
   VBI(3)=MXP
   VBI(4)=1  ;***0 = Upper tail as in Drezner & Wesolowsky; 1 = Bottom tail***;
   VBI(5)=1  ;***0 = 3 pt approximation; 1 = 5 point approximation***;
   BV=BIVARIATE(VBI)
   JP=((DV-1)*(DVP-1)+(DV-1)*(1-2*DVP)*PHIMXP+(DVP-1)*(1-2*DV)*PHIMX+(1-2*DV)*(1-2*DVP)*BV)
  ENDIF
   IF(JP.LE.0.0) EXIT
   LOGL=LOG(JP/PCP)
   Y = -2*LOGL

  MXP=MX
  PCP=PC
  DVP=DV
  TIMEP=TIME
  PHIMXP=PHIMX

$THETA     
      -1.7   ; 1  B1
       1.2   ; 2  B2
       2.9   ; 3  B3
       1.4   ; 4  LOG(B4)
       1.2   ; 5  LOG(B5)
       (0.0 FIXED)   ; 6  LOG SQRT VAR(ETA1)
       (0.0 FIXED)  ; 7  LOG SQRT VAR(ETA2)
       2.2   ; 8  RHO parameter

$OMEGA DIAGONAL(2)
      0.8       ; V1
      0.8       ; V2

;$EST METHOD=IMP LAPLACE -2LL PRINT=1 NITER=300 ISAMPLE=300 SIGL=6 CTYPE=3 NOHABORT
$EST MAX=0 PRINT=1 METHOD=1 LAPLACE -2LL SIGL=10 NOHABORT
$COV COMPRESS MATRIX=R PRINT=E UNCONDITIONAL