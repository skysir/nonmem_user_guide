Mon 01/27/2020 
04:24 PM
;Model Desc: Receptor Mediated Clearance model with Dynamic Change in Receptors
;Project Name: antibody
;Project ID: NO PROJECT DESCRIPTION

;Project ID: NO PROJECT DESCRIPTION

$PROB RUN# wexample6x (from r2compl), using FAST
$ABBR DERIV2=NO DERIV2=NOCOMMON ; DERIV1=NO
$INPUT C SET ID JID TIME DV AMT RATE EVID MDV CMT
$DATA a_uflg.csv IGNORE=C IGNORE=(AMT.EQN.1.0) IGNORE=(RATE.EQN.1.0E+08)

; The new numerical integration solver is used, although ADVAN=9 is also efficient
; for this problem.
$SUBROUTINES ADVAN13 TRANS1 TOL=6 ATOL=6
$MODEL NCOMPARTMENTS=3

;Initial Omegas
$OMEGA BLOCK(8) VALUES(0.5,0.01)

$PK
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
MU_5=THETA(5)
MU_6=THETA(6)
MU_7=THETA(7)
MU_8=THETA(8)
VC=EXP(MU_1+ETA(1))
K10=EXP(MU_2+ETA(2))
K12=EXP(MU_3+ETA(3))
K21=EXP(MU_4+ETA(4))
VM=EXP(MU_5+ETA(5))
KMC=EXP(MU_6+ETA(6))
K03=EXP(MU_7+ETA(7))
K30=EXP(MU_8+ETA(8))
S3=VC
S1=VC
KM=KMC*S1
A_0(3)=K03/K30
MTIME(1)=7.1
MTDIFF=1
AZTEST=A_0FLG
IF(TSTATE==MTIME(1).AND.AZTEST==0) A_UFLG=1
A_U(1)=A(1)+1000.0
A_U(2)=A(2)
A_U(3)=A(3)

$DES
DADT(1) = -(K10+K12)*A(1) + K21*A(2) - VM*A(1)*A(3)/(A(1)+KM)
DADT(2) = K12*A(1) - K21*A(2)
DADT(3) =  -(VM-K30)*A(1)*A(3)/(A(1)+KM) - K30*A(3) + K03

$ERROR
ETYPE=1
IF(CMT.NE.1) ETYPE=0
CP=A(1)/S1
CR=A(3)/S3
IPRE=CP
IF(CMT.NE.1) IPRE=CR
Y = IPRE + IPRE*ETYPE*EPS(1) + IPRE*(1.0-ETYPE)*EPS(2)


$THETA 
;Initial Thetas
( 4.0 )  ;[MU_1]
( -3.1 ) ;[MU_2]
( 0.5 )  ;[MU_3]
( -0.2 );[MU_4]      
( 3.2 ) ;[MU_5]
( 0.01 )  ;[MU_6]
( 4.0 )  ;[MU_7]
( -0.1) ;[MU_8]



$SIGMA  
0.1 ;[p]
0.1 ;[p]

$EST METHOD=1 INTERACTION PRINT=1 NOABORT NOPRIOR=1 SIGL=6 MAXEVAL=0 NSIG=2 MCETA=10 
     NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST FAST
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       27 JAN 2020
Days until program expires :3775
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0 Beta version 3
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 RUN# wexample6x (from r2compl), using FAST
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     1650
 NO. OF DATA ITEMS IN DATA SET:  11
 ID DATA ITEM IS DATA ITEM NO.:   3
 DEP VARIABLE IS DATA ITEM NO.:   6
 MDV DATA ITEM IS DATA ITEM NO.: 10
0INDICES PASSED TO SUBROUTINE PRED:
   9   5   7   8   0   0  11   0   0   0   0
0LABELS FOR DATA ITEMS:
 C SET ID JID TIME DV AMT RATE EVID MDV CMT
0FORMAT FOR DATA:
 (2E2.0,2E3.0,E5.0,E10.0,E5.0,4E2.0)

 TOT. NO. OF OBS RECS:     1568
 TOT. NO. OF INDIVIDUALS:       50
0LENGTH OF THETA:   8
0DEFAULT THETA BOUNDARY TEST OMITTED:   YES
0OMEGA HAS BLOCK FORM:
  1
  1  1
  1  1  1
  1  1  1  1
  1  1  1  1  1
  1  1  1  1  1  1
  1  1  1  1  1  1  1
  1  1  1  1  1  1  1  1
0DEFAULT OMEGA BOUNDARY TEST OMITTED:   YES
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   2
0DEFAULT SIGMA BOUNDARY TEST OMITTED:   YES
0INITIAL ESTIMATE OF THETA:
   0.4000E+01 -0.3100E+01  0.5000E+00 -0.2000E+00  0.3200E+01  0.1000E-01  0.4000E+01 -0.1000E+00
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.5000E+00
                  0.1000E-01   0.5000E+00
                  0.1000E-01   0.1000E-01   0.5000E+00
                  0.1000E-01   0.1000E-01   0.1000E-01   0.5000E+00
                  0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.5000E+00
                  0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.5000E+00
                  0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.5000E+00
                  0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.1000E-01   0.5000E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+00
 0.0000E+00   0.1000E+00
1DOUBLE PRECISION PREDPP VERSION 7.5.0 Beta version 3

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF EQUATIONS (LSODA, ADVAN13)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   7
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         COMP 1       ON         YES        YES        YES        YES
    2         COMP 2       ON         YES        YES        NO         NO
    3         COMP 3       ON         YES        YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
 INITIAL (BASE) TOLERANCE SETTINGS:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:   6
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            9           *           *           *           *
    2            *           *           *           *           *
    3            8           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0FIRST MODEL TIME PARAMETER ASSIGNED TO ROW NO.: 10
 LAST  MODEL TIME PARAMETER ASSIGNED TO ROW NO.: 10
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      9
   TIME DATA ITEM IS DATA ITEM NO.:          5
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   7
   DOSE RATE DATA ITEM IS DATA ITEM NO.:     8
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:   11

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE TIMES AND AT MODEL TIMES.

0PK SUBROUTINE INDICATES THAT COMPARTMENT AMOUNTS ARE INITIALIZED.
0PK SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 
 
 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction (No Prior) (Evaluation)
 
 ESTIMATION STEP OMITTED:                 YES
 ANALYSIS TYPE:                           POPULATION
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    10
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      6
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     6
 NOPRIOR SETTING (NOPRIOR):                 1
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): a_uflg.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:   6
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:   6
 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 Elapsed evaluation time in seconds:     2.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************   FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR) (EVALUATION)  ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3110.228       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************   FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR) (EVALUATION)  ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8     
 
         4.00E+00 -3.10E+00  5.00E-01 -2.00E-01  3.20E+00  1.00E-02  4.00E+00 -1.00E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        5.00E-01
 
 ETA2
+        1.00E-02  5.00E-01
 
 ETA3
+        1.00E-02  1.00E-02  5.00E-01
 
 ETA4
+        1.00E-02  1.00E-02  1.00E-02  5.00E-01
 
 ETA5
+        1.00E-02  1.00E-02  1.00E-02  1.00E-02  5.00E-01
 
 ETA6
+        1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.00E-02  5.00E-01
 
 ETA7
+        1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.00E-02  5.00E-01
 
 ETA8
+        1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.00E-02  5.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1      EPS2     
 
 EPS1
+        1.00E-01
 
 EPS2
+        0.00E+00  1.00E-01
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        7.07E-01
 
 ETA2
+        2.00E-02  7.07E-01
 
 ETA3
+        2.00E-02  2.00E-02  7.07E-01
 
 ETA4
+        2.00E-02  2.00E-02  2.00E-02  7.07E-01
 
 ETA5
+        2.00E-02  2.00E-02  2.00E-02  2.00E-02  7.07E-01
 
 ETA6
+        2.00E-02  2.00E-02  2.00E-02  2.00E-02  2.00E-02  7.07E-01
 
 ETA7
+        2.00E-02  2.00E-02  2.00E-02  2.00E-02  2.00E-02  2.00E-02  7.07E-01
 
 ETA8
+        2.00E-02  2.00E-02  2.00E-02  2.00E-02  2.00E-02  2.00E-02  2.00E-02  7.07E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1      EPS2     
 
 EPS1
+        3.16E-01
 
 EPS2
+        0.00E+00  3.16E-01
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,        3.744
Stop Time: 
Mon 01/27/2020 
04:24 PM
