;Model Desc:  ar1est_01mod2.ctl
;Project Name: testprob.08.08.08
;Project ID: None

$PROB  ctlar1mod same as ar1est_01mod2.ctl uses new features
$ABBR declare T(NO)
$ABBR DECLARE DOWHILE J
$ABBR declare integer i
$INPUT CX ID DOSE=AMT TIME CP=DV WT MDV
$DATA       ar1sim01mod.dat IGNORE=@

$SUBROUTINES  ADVAN2

$PK
;THETA(1)=MEAN ABSORPTION RATE CONSTANT (1/HR)
;THETA(2)=MEAN ELIMINATION RATE CONSTANT (1/HR)
;THETA(3)=SLOPE OF CLEARANCE VS WEIGHT RELATIONSHIP (LITERS/HR/KG)
;SCALING PARAMETER=VOLUME/WT SINCE DOSE IS WEIGHT-ADJUSTED
   CALLFL=1
   KA=THETA(1)*DEXP(ETA(1))
   K=THETA(2)*DEXP(ETA(2))
   CL=THETA(3)*WT*DEXP(ETA(3))
   SC=CL/K/WT

$ERROR
 IF(NEWIND.NE.2) I=0
 IF(MDV.EQ.0)THEN
    I=I+1
    T(I)=TIME
    J=1
    DO WHILE (J<=I)
    CORRL2(J,1)=EXP(-THETA(4)*(TIME-T(J)))
    J=J+1
    ENDDO
 ENDIF
   Y=F+F*EPS(1)

$THETA  
(.1,3,5)      ;[KA]
(.008,.08,.5) ;[KE]
(.004,.04,.9) ;[CL]
(0,.5,1)      ;[Ar_parameter]


$OMEGA
 .09      ;[P]
$OMEGA BLOCK(2)
 .09      ;[P]
 .04      ;[F]
 .09      ;[P]

$SIGMA  
 .2       ;[P]

$EST  
  NSIG=3 
  METHOD=1
  INTERACTION
  MAXEVAL=9999
  PRINT=5
  NOTHETABOUNDTEST
  NOOMEGABOUNDTEST
  NOSIGMABOUNDTEST 
  MSFO=ar1est_01mod2.msf

$COV MATRIX=S PRINT=E

$TABLE ID DOSE TIME CP WT 
  NOPRINT
  ONEHEADER
  FILE=tabar1mod

