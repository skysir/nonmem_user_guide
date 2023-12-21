Thu Jan 22 08:56:32 PST 2015
$PROBLEM  Parent drug, using ADVAN6 
$INPUT ID TIME AMT UVOL DV CMT MDV EVID L2
$DATA dextroparent.dat IGNORE=#
$SUBROUTINES ADVAN6 TRANS1 TOL=4
$MODEL
  COMP=(DEPOT)
  COMP=(PLASMA DEFOBS) ; PARENT IN PLASMA
$PK
  K12=THETA(1)
  MU1=LOG(THETA(2))
  V2=EXP(MU1+ETA(1))
  MU2=LOG(THETA(3))
  CLP=EXP(MU2+ETA(2)) ; RENAL CL FOR PARENT
  CLB=THETA(4)         ; METABOLIC CL FOR METABOLIC
  K23=CLP/V2+CLB/V2
  F0=CLP/(CLP+CLB)
  S2=V2
  S3=UVOL
$ERROR (EVERY EVENT)
  ACMT=ABS(CMT) ; output compartment may have negative value of CMT
  IF(ACMT.EQ.2) Y=F*(1+EPS(1))
  IF(ACMT.EQ.3) Y=F*(1+EPS(2))
$DES
  DADT(1)=-K12*A(1)
  DADT(2)=K12*A(1)-K23*A(2)
$THETA
  (0.01,0.8,6) ;KA
  (0.01,43,1000);V2
  (0.0001,20,190);CLP
  (0.01,15,90);CLB
$OMEGA
  0.05  0.05
$SIGMA 
  .01 .01
  .01 .01 ; eps(3) and eps(4) for consistent EPS with sim_metab_a6 example
$SIM (111111) ONLYSIM
$TABLE ID TIME AMT UVOL DV SIMP=PRED CMT MDV EVID L2 NOAPPEND FILE=sim_parent_a66.tab
 
 
NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  3) THERE MAY BE AN ERROR IN THE ABBREVIATED CODE. THE FOLLOWING
 ONE OR MORE RANDOM VARIABLES ARE DEFINED WITH "IF" STATEMENTS THAT DO NOT
 PROVIDE DEFINITIONS FOR BOTH THE "THEN" AND "ELSE" CASES. IF ALL
 CONDITIONS FAIL, THE VALUES OF THESE VARIABLES WILL BE ZERO.
  
   Y

             
 (WARNING  32) $SIGMA INCLUDES A NON-FIXED INITIAL ESTIMATE CORRESPONDING TO
 AN EPS (OR ERR) THAT IS NOT USED IN ABBREVIATED CODE.
 
 
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM)    DOUBLE PRECISION NONMEM    VERSION VI LEVEL 2.0  
 DEVELOPED AND PROGRAMMED BY STUART BEAL AND LEWIS SHEINER
 
 PROBLEM NO.:         1
 Parent drug, using ADVAN6                                               
0DATA CHECKOUT RUN:              NO 
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO 
 NO. OF DATA RECS IN DATA SET:        7
 NO. OF DATA ITEMS IN DATA SET:   9
 ID DATA ITEM IS DATA ITEM NO.:   1
 L2 DATA ITEM IS DATA ITEM NO.:   9
 DEP VARIABLE IS DATA ITEM NO.:   5
 MDV DATA ITEM IS DATA ITEM NO.:  7
0INDICES PASSED TO SUBROUTINE PRED:
  8  2  3  0  0  0  6  0  0
  0  0
0LABELS FOR DATA ITEMS:
   ID    TIME     AMT    UVOL      DV     CMT     MDV    EVID      L2
0LABELS FOR SPECIAL ITEMS:
 SIMP     RES    WRES
0FORMAT FOR DATA:
 (2E6.0,E10.0,E8.0,E6.0,E7.0,2E6.0,E2.0)                                         
 
 TOT. NO. OF OBS RECS:        4
 TOT. NO. OF INDIVIDUALS:      1
0LENGTH OF THETA:  4
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO 
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:  2
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO 
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:  4
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO 
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.1000E-01     0.8000E+00     0.6000E+01
  0.1000E-01     0.4300E+02     0.1000E+04
  0.1000E-03     0.2000E+02     0.1900E+03
  0.1000E-01     0.1500E+02     0.9000E+02
0INITIAL ESTIMATE OF OMEGA:
 0.5000E-01
 0.0000E+00   0.5000E-01
0INITIAL ESTIMATE OF SIGMA:
 0.1000E-01
 0.0000E+00   0.1000E-01
 0.0000E+00   0.0000E+00   0.1000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.1000E-01
0SIMULATION STEP OMITTED:    NO 
 OBJ FUNC EVALUATED:         NO 
 SOURCE  1:
    SEED1:        111111   SEED2:             0   PSEUDO-NORMAL       
0TABLES STEP OMITTED:    NO 
 NO. OF TABLES:           1
0-- TABLE  1 --
04 COLUMNS APPENDED:     NO 
 PRINTED:               YES 
 FOR TABLE FILE,
 HEADER:                YES 
 FILE TO BE FORWARDED:   NO 
0USER-CHOSEN ITEMS 
 IN THE ORDER THEY WILL APPEAR IN THE TABLE:
   ID    TIME     AMT    UVOL      DV    SIMP     CMT     MDV    EVID      L2
1DOUBLE PRECISION PREDPP VERSION V LEVEL 2.0  
 
 GENERAL NONLINEAR KINETICS MODEL (ADVAN6)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   2
0COMPARTMENT ATTRIBUTES 
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        ON         YES        YES        YES        NO 
    2         PLASMA       ON         YES        YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO 
0NRD VALUE FROM SUBROUTINE TOL:   4
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG  
    1           *           *           *           *           *
    2           4           *           *           *           *
    3           5           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0OUTPUT FRACTION PARAMETER ASSIGNED TO ROW NO.:  3
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      8
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   3
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:    6
 
0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 SIMULATION STEP PERFORMED
 SOURCE  1:
    SEED1:    1422409161   SEED2:             0
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          TABLES OF DATA AND PREDICTIONS                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 
1TABLE NO.  1



 LINE NO.     ID      TIME       AMT      UVOL        DV      SIMP       CMT       MDV      EVID        L2
 
    1
+        2.00E+00  0.00E+00  3.00E+04  0.00E+00  0.00E+00  0.00E+00  1.00E+00  1.00E+00  1.00E+00  1.00E+00
 
    2
+        2.00E+00  1.71E-01  0.00E+00  0.00E+00  0.00E+00  8.31E+01  3.00E+00  1.00E+00  2.00E+00  1.00E+00
 
    3
+        2.00E+00  2.00E+00  0.00E+00  0.00E+00  2.43E+02  2.22E+02  2.00E+00  0.00E+00  0.00E+00  1.00E+00
 
    4
+        2.00E+00  2.00E+00  0.00E+00  9.31E+01  9.61E+01  8.67E+01 -3.00E+00  0.00E+00  0.00E+00  1.00E+00
 
    5
+        2.00E+00  2.00E+00  0.00E+00  0.00E+00  0.00E+00  2.22E+02  3.00E+00  1.00E+00  2.00E+00  2.00E+00
 
    6
+        2.00E+00  3.15E+00  0.00E+00  0.00E+00  1.68E+02  1.38E+02  2.00E+00  0.00E+00  0.00E+00  2.00E+00
 
    7
+        2.00E+00  3.15E+00  0.00E+00  1.34E+02  2.75E+01  3.10E+01 -3.00E+00  0.00E+00  0.00E+00  2.00E+00
 
Thu Jan 22 08:56:32 PST 2015