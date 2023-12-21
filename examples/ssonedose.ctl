; SSONEDOSE.CTL 
;
$PROBLEM DOSE SUPERPOSITION WITH STEADY STATE, USING AN SS RECORD.
; THIS IS AN EXAMPLE OF STEADY STATE WITH MULTIPLE BOLUS  DOSES 
; INTO A TRANSIT COMPARTMENT.
; COMPARE OUTPUT WITH SSMULTIDOSE.CTL
$INPUT ID TIME CMT AMT DV II SS
$DATA ssonedose.dat IGNORE=@
$SUBROUTINES ADVAN6 TOL=10 TRANS1

$MODEL COMP=ABS; (DEFDOSE)
       COMP=CENTRAL ;(DEFOBS)
       COMP=PERI
$ABBR DECLARE DOWHILE I
$ABBR DECLARE REAL DOSETIME,DOSE,DELTAT,ABI,DESTOL

$PK 
CALLFL=-2
DESTOL=10**(-17) ; DES SHOULD COMPUTE MORE SIG. DIGITS THAN THE TOL VALUE

IF (SS == 1 ) THEN  ; SS DOSE RECORD
DOSETIME=TIME
DOSE=AMT
IIVAL=II
ENDIF

CL = THETA(1)*EXP(ETA(1))
V2 = THETA(2)*EXP(ETA(2))
Q  = THETA(3)
V3 = THETA(4)

K=CL/V2
S2 =V2
K23=Q/V2
K32=Q/V3

;--------ABSORPTION MODEL-----------
F1=0
MTT  = THETA(5)*EXP(ETA(4))
NN   = THETA(6)*EXP(ETA(5))
BIO  = THETA(7)*EXP(ETA(3))
KTR  = (NN+1)/MTT

NFAC=SQRT(2*3.1416)*NN**(NN+0.5)*(EXP(-NN))*(1+1/(12*NN))
KINPT=BIO*KTR**(NN+1)/NFAC

$DES
; DES ADDS DRUG FROM IMPLIED DOSES GOING BACK IN TIME TILL
; THE AMOUNT IS NEGLIGIBLE (ABS(IPT) < DESTOL)
; THE VALUE OF DESTOL WAS SET EXPERIMENTALLY, TO AVOID ERROR MESSAGE
; NUMERICAL DIFFICULTIES WITH INTEGRATION ROUTINE
; FROM PREDPP.
; THE COUNTER I IS USED TO PREVENT A RUNAWAY SITUATION.
ABI=1
INPT=0
DELTAT=T-DOSETIME-IIVAL
I=1
DOWHILE (I < 999 .AND. ABI >= DESTOL )
IPT=0
; IT IS IMPORTANT TO COMPUTE DELTAT HERE, BEFORE THE IF/ENDIF
DELTAT=DELTAT+IIVAL
IF (DELTAT>0) THEN
IPT=DOSE*DELTAT**NN*EXP(-KTR*DELTAT)
INPT=INPT+IPT
ABI=ABS(IPT)
ENDIF
I=I+1
ENDDO

 DADT(1)=KINPT*INPT-KTR*A(1)
 DADT(2)=KTR*A(1)-K23*A(2)-K*A(2)+K32*A(3)
 DADT(3)=K23*A(2)-K32*A(3)

$ERROR

Y = F*EXP(EPS(1))

$THETA 
4 ;CL
3 ;V2
1 ;Q
2 ;V3
10 ;MTT
5 ;NN
0.85 ;BIO
$OMEGA 1 1 1 1 1 
$SIGMA 1 
$TABLE  ID TIME FORMAT=S1PE18.11 FILE=ssonedose.tab
