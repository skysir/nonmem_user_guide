Mon 09/30/2013 
03:05 PM
; Modified from appendix 4 of Jun Shen, Alison Boeckmann & Andrew Vick
; Modified for NONMEM 7.3
; Journal of Pharmacokinetics and Pharmacodynamics
; ISSN 1567-567X Volume 39 Number 3
; J Pharmacokinet Pharmacodyn (2012) 39:251-262 DOI 10.1007/s10928-012-9247-3

$PROBLEM Dose Superposition
$INPUT ID TIME CMT AMT DV ADDL II RATE
$DATA sumdosetf.csv IGNORE=@  ; identical to sall6.csv
$SUBROUTINES ADVAN6 TOL=10 TRANS1

$MODEL COMP=ABS ;(DEFDOSE)
       COMP=CENTRAL;(DEFOBS)
       COMP=PERI

$abbr declare dosetime(100),dose(100)
$abbr declare dowhile i
$abbr declare dowhile ndose

$PK 
CALLFL=-2
IF (NEWIND < 2) NDOSE=0

IF (AMT > 0 .and. cmt==1) THEN
 NDOSE=NDOSE+1
 dosetime(NDOSE)=TIME
 DOSE(NDOSE)=AMT
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
INPT=0
I=1
DOWHILE (I<=NDOSE)
IPT=0
IF (T>=dosetime(I)) IPT=DOSE(I)*(T-dosetime(I))**NN*EXP(-KTR*(T-dosetime(I)))
INPT=INPT+IPT
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

$OMEGA 0 FIX 0 FIX 0 FIX 0 FIX 0 FIX 
$SIGMA 0 FIX
$SIMULATION (123456) ONLY
;$ESTIMATION METHOD=0
$TABLE  ID TIME CMT CL V2 Q V3 MTT NN INPT ONEHEADER NOPRINT FORMAT=s1PE18.11 FILE=sumdosetn.txt
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  98) A RANDOM QUANTITY IS RAISED TO A POWER. IF THE RESULT AFFECTS
 THE VALUE OF THE OBJECTIVE FUNCTION AND THE "ONLYSIM" OPTION IS NOT USED,
 THE USER SHOULD ENSURE THAT THE RANDOM QUANTITY IS NEVER 0 WHEN THE POWER
 IS < 1.
             
 (WARNING  100) A RANDOM QUANTITY IS USED AS A POWER. IF THE RESULT AFFECTS
 THE VALUE OF THE OBJECTIVE FUNCTION AND THE "ONLYSIM" OPTION IS NOT USED,
 THE USER SHOULD ENSURE THAT THE QUANTITY RAISED TO THE POWER IS NOT 0.

 (DATA WARNING   2) RECORD         1, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         2, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         3, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         4, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         5, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         6, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         7, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         8, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD         9, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        10, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        11, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        12, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        13, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        14, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        15, DATA ITEM   8, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        16, DATA ITEM   4, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        17, DATA ITEM   4, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        18, DATA ITEM   4, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        19, DATA ITEM   4, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.

 (DATA WARNING   2) RECORD        20, DATA ITEM   5, CONTENTS:  
 THE NUMBER OF DATA ITEMS SPECIFIED IN $INPUT EXCEEDS THE NUMBER OF VALUES
 IN A RECORD OF THE NM-TRAN DATA FILE. NULLS WERE SUPPLIED FOR MISSING
 VALUES, STARTING WITH THE ABOVE NUMBERED DATA ITEM.*
             
 (WARNING  48) DES-DEFINED ITEMS ARE COMPUTED ONLY WHEN EVENT TIME
 INCREASES. E.G., DISPLAYED VALUES ASSOCIATED WITH THE FIRST EVENT RECORD
 OF AN INDIVIDUAL RECORD ARE COMPUTED WITH (THE LAST ADVANCE TO) AN EVENT
 TIME OF THE PRIOR INDIVIDUAL RECORD.

 * THE MAXIMUM NUMBER OF WARNINGS OF ONE OR MORE TYPES WAS REACHED.
 IT IS POSSIBLE THAT SOME WARNING MESSAGES WERE SUPPRESSED.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       30 SEP 2013
Days until program expires :6087
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.3beta7.0(P)
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 Dose Superposition
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      125
 NO. OF DATA ITEMS IN DATA SET:  10
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   5
 MDV DATA ITEM IS DATA ITEM NO.: 10
0INDICES PASSED TO SUBROUTINE PRED:
   9   2   4   8   0   7   3   0   0   0   6
0LABELS FOR DATA ITEMS:
 ID TIME CMT AMT DV ADDL II RATE EVID MDV
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 CL V2 Q V3 MTT NN INPT
0FORMAT FOR DATA:
 (8E6.0,2F2.0)

 TOT. NO. OF OBS RECS:      119
 TOT. NO. OF INDIVIDUALS:      1
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
   0.4000E+01  0.3000E+01  0.1000E+01  0.2000E+01  0.1000E+02  0.5000E+01  0.8500E+00
0INITIAL ESTIMATE OF OMEGA:
 0.0000E+00
 0.0000E+00   0.0000E+00
 0.0000E+00   0.0000E+00   0.0000E+00
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.0000E+00
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0SIMULATION STEP OMITTED:    NO
 OBJ FUNC EVALUATED:         NO
0SIMULATION RANDOM METHOD SELECTED (RANMETHOD): DEFAULT
 SOURCE   1:
   SEED1:        123456   SEED2:             0   PSEUDO-NORMAL
0WARNING: NO. OF OBS RECS IN INDIVIDUAL REC NO.      1 (IN INDIVIDUAL
 REC ORDERING) EXCEEDS ONE WHILE INITIAL ESTIMATE OF WITHIN INDIVIDUAL VARIANCE IS ZERO
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    11456
 RANMETHOD:
 MC SAMPLES (ESEED):    300
 WRES SQUARE ROOT TYPE:            EIGENVALUE
0-- TABLE   1 --
 PRINTED:                NO
 HEADER:                YES
 FILE TO BE FORWARDED:   NO
 FORMAT:                s1PE18.11
 LFORMAT:
 RFORMAT:
0USER-CHOSEN ITEMS:
 ID TIME CMT CL V2 Q V3 MTT NN INPT
1DOUBLE PRECISION PREDPP VERSION 7.3beta7.0(P)

 GENERAL NONLINEAR KINETICS MODEL (ADVAN6)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   6
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         ABS          ON         YES        YES        YES        NO
    2         CENTRAL      ON         YES        YES        NO         YES
    3         PERI         ON         YES        YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
0NRD VALUE FROM SUBROUTINE TOL:  10
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           8           *           *           *
    2            7           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0ERROR IN LOG Y IS MODELED
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      9
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4
   DOSE RATE DATA ITEM IS DATA ITEM NO.:     8
   INTERVAL DATA ITEM IS DATA ITEM NO.:      7
   ADDL. DOSES DATA ITEM IS DATA ITEM NO.:   6
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:    3

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE TIMES.
0DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 OTHERWISE, ERROR SUBROUTINE CALLED ONCE IN THIS PROBLEM.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 SIMULATION STEP PERFORMED
 SOURCE  1:
    SEED1:     345794978   SEED2:             0
 ESTIMATION STEP OMITTED:                 YES 
 #CPUT: Total CPU Time in Seconds,        0.125
Stop Time: 
Mon 09/30/2013 
03:05 PM
