Mon 10/31/2016 
12:18 AM
;Model Desc: Receptor Mediated Clearance model with Dynamic Change in Receptors
;Project Name: nm7examples
;Project ID: NO PROJECT DESCRIPTION

$PROB RUN# example6 (from r2compl)
$INPUT C SET ID JID TIME DV=CONC DOSE=AMT RATE EVID MDV CMT
$DATA example6.csv IGNORE=C

; The new numerical integration solver is used, although ADVAN=9 is also efficient
; for this problem.
$SUBROUTINES ADVAN13 TRANS1 TOL=4
$MODEL NCOMPARTMENTS=3


$PK
include nonmem_reserved_general
MUFIRSTREC=1
OBJQUICK=2
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
F3=K03/K30

$DES
DADT(1) = -(K10+K12)*A(1) + K21*A(2) - VM*A(1)*A(3)/(A(1)+KM)
DADT(2) = K12*A(1) - K21*A(2)
DADT(3) =  -VM*A(1)*A(3)/(A(1)+KM) - K30*A(3) + K03

$ERROR
CALLFL=0
ETYPE=1
IF(CMT.NE.1) ETYPE=0
IPRED=F
Y = F + F*ETYPE*EPS(1) + F*(1.0-ETYPE)*EPS(2)


;Initial Thetas
$THETA
( 4.0 )  ;[MU_1]
( -2 ) ;[MU_2]
( 1.0 )  ;[MU_3]
( -0.2 );[MU_4]      
( 2.2 ) ;[MU_5]
( 0.5 )  ;[MU_6]
( 4.0 )  ;[MU_7]
( -0.8) ;[MU_8]

;Initial Omegas
$OMEGA BLOCK(8) VALUES(0.8,0.001)

$SIGMA  
0.1 ;[p]
0.1 ;[p]

$PRIOR NWPRI
$OMEGAP BLOCK(8) FIXED VALUES(0.2,0.0)
$OMEGAPD 8.0 FIXED

$EST METHOD=ITS INTERACTION SIGL=4 NITER=15 PRINT=1 NOABORT NOPRIOR=1 file=example6hmto19_its.ext
$EST METHOD=bayes INTERACTION NBURN=2000 NITER=0 PRINT=10 MASSRESET=1 NOPRIOR=0 file=example6hmto19_bayes.ext
$EST METHOD=NUTS INTERACTION  NBURN=250 NITER=2000 PRINT=5 MASSRESET=0 PMADAPT=200  file=example6hmto19.ext
     OLKJDF=8.0
$COV MATRIX=R UNCONDITIONAL
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

 (MU_WARNING 20) MU_001: MU_ VARIABLE SHOULD NOT BE DEFINED AFTER VERBATIM CODE.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       31 OCT 2016
Days until program expires :4961
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.4.0 alpha12 (nm74a12)
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 RUN# example6 (from r2compl)
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     1750
 NO. OF DATA ITEMS IN DATA SET:  11
 ID DATA ITEM IS DATA ITEM NO.:   3
 DEP VARIABLE IS DATA ITEM NO.:   6
 MDV DATA ITEM IS DATA ITEM NO.: 10
0INDICES PASSED TO SUBROUTINE PRED:
   9   5   7   8   0   0  11   0   0   0   0
0LABELS FOR DATA ITEMS:
 C SET ID JID TIME CONC DOSE RATE EVID MDV CMT
0FORMAT FOR DATA:
 (2E2.0,2E3.0,E5.0,E10.0,2E5.0,3E2.0)

 TOT. NO. OF OBS RECS:     1568
 TOT. NO. OF INDIVIDUALS:     50
0LENGTH OF THETA:   9
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS BLOCK FORM:
  1
  1  1
  1  1  1
  1  1  1  1
  1  1  1  1  1
  1  1  1  1  1  1
  1  1  1  1  1  1  1
  1  1  1  1  1  1  1  1
  0  0  0  0  0  0  0  0  2
  0  0  0  0  0  0  0  0  2  2
  0  0  0  0  0  0  0  0  2  2  2
  0  0  0  0  0  0  0  0  2  2  2  2
  0  0  0  0  0  0  0  0  2  2  2  2  2
  0  0  0  0  0  0  0  0  2  2  2  2  2  2
  0  0  0  0  0  0  0  0  2  2  2  2  2  2  2
  0  0  0  0  0  0  0  0  2  2  2  2  2  2  2  2
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   2
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
 -0.1000E+07     0.4000E+01     0.1000E+07
 -0.1000E+07    -0.2000E+01     0.1000E+07
 -0.1000E+07     0.1000E+01     0.1000E+07
 -0.1000E+07    -0.2000E+00     0.1000E+07
 -0.1000E+07     0.2200E+01     0.1000E+07
 -0.1000E+07     0.5000E+00     0.1000E+07
 -0.1000E+07     0.4000E+01     0.1000E+07
 -0.1000E+07    -0.8000E+00     0.1000E+07
  0.8000E+01     0.8000E+01     0.8000E+01
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.8000E+00
                  0.1000E-02   0.8000E+00
                  0.1000E-02   0.1000E-02   0.8000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.8000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.8000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.8000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.8000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.8000E+00
        2                                                                                  YES
                  0.2000E+00
                  0.0000E+00   0.2000E+00
                  0.0000E+00   0.0000E+00   0.2000E+00
                  0.0000E+00   0.0000E+00   0.0000E+00   0.2000E+00
                  0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.2000E+00
                  0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.2000E+00
                  0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.2000E+00
                  0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.2000E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+00
 0.0000E+00   0.1000E+00
0COVARIANCE STEP OMITTED:        NO
 R MATRIX SUBSTITUTED:          YES
 S MATRIX SUBSTITUTED:           NO
 EIGENVLS. PRINTED:              NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 TURN OFF Cholesky Transposition of R Matrix (CHOLROFF): NO
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):              -1
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 0
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING: (FPOSDEF):0
0
 PRIOR SUBROUTINE USER-SUPPLIED
1DOUBLE PRECISION PREDPP VERSION 7.4.0 alpha12 (nm74a12)

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF EQUATIONS (ADVAN13)
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
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            9           *           *           *           *
    2            *           *           *           *           *
    3            8          10           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      9
   TIME DATA ITEM IS DATA ITEM NO.:          5
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   7
   DOSE RATE DATA ITEM IS DATA ITEM NO.:     8
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:   11

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 OTHERWISE, ERROR SUBROUTINE CALLED ONLY WITH OBSERVATION EVENTS.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 
 
 #TBLN:      1
 #METH: Iterative Two Stage (No Prior)
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            3480
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      4
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     4
 NOPRIOR SETTING (NOPRIOR):                 ON
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example6hmto19_its.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 WISHART PRIOR DF INTERPRETATION (WISHTYPE):0
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 EM OR BAYESIAN METHOD USED:                ITERATIVE TWO STAGE (ITS)
 MU MODELING PATTERN (MUM):
 GRADIENT/GIBBS PATTERN (GRD):
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  0
 ITERATIONS (NITER):                        15
 ANEAL SETTING (CONSTRAIN):                 1

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8
 THETAS THAT ARE SIGMA-LIKE:
 
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=  -3179.34804709958
 iteration            1 OBJ=  -3563.39825052570
 iteration            2 OBJ=  -3698.49131024100
 iteration            3 OBJ=  -3813.49052442310
 iteration            4 OBJ=  -3921.38645757702
 iteration            5 OBJ=  -4025.66889419349
 iteration            6 OBJ=  -4127.61073505300
 iteration            7 OBJ=  -4227.50665677951
 iteration            8 OBJ=  -4325.17010194533
 iteration            9 OBJ=  -4419.82458046677
 iteration           10 OBJ=  -4509.87267715267
 iteration           11 OBJ=  -4592.19272921940
 iteration           12 OBJ=  -4659.87003784235
 iteration           13 OBJ=  -4700.24430332201
 iteration           14 OBJ=  -4709.62334419566
 iteration           15 OBJ=  -4710.51279864924
 
 #TERM:
 OPTIMIZATION WAS NOT TESTED FOR CONVERGENCE


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -8.7799E-04 -1.6780E-03  2.3470E-03  1.3344E-03  1.1374E-03  1.5436E-03  2.9935E-04  8.2864E-04
 SE:             6.9351E-02  5.3938E-02  3.7889E-02  6.5123E-02  5.6761E-02  5.7438E-02  6.4405E-02  6.1401E-02
 N:                      50          50          50          50          50          50          50          50
 
 P VAL.:         9.8990E-01  9.7518E-01  9.5061E-01  9.8365E-01  9.8401E-01  9.7856E-01  9.9629E-01  9.8923E-01
 
 ETASHRINKSD(%)  6.4398E-01  4.5093E+00  8.2195E+00  1.6470E+00  1.4992E+00  5.8680E+00  3.9429E-01  1.5720E+00
 ETASHRINKVR(%)  1.2838E+00  8.8152E+00  1.5763E+01  3.2670E+00  2.9760E+00  1.1392E+01  7.8702E-01  3.1193E+00
 EBVSHRINKSD(%)  6.3658E-01  5.3912E+00  9.8831E+00  2.1078E+00  1.5741E+00  6.3487E+00  4.5677E-01  1.7624E+00
 EBVSHRINKVR(%)  1.2691E+00  1.0492E+01  1.8789E+01  4.1712E+00  3.1235E+00  1.2294E+01  9.1146E-01  3.4938E+00
 EPSSHRINKSD(%)  1.5721E+01  7.2546E+00
 EPSSHRINKVR(%)  2.8970E+01  1.3983E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         1568
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    2881.79124012985     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -4710.51279864924     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1828.72155851938     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           400
  
 #TERE:
 Elapsed estimation  time in seconds:    46.32
 Elapsed covariance  time in seconds:     0.31
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -4710.513       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8     
 
         3.91E+00 -2.20E+00  5.57E-01 -1.86E-01  2.26E+00  2.14E-01  3.71E+00 -7.07E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        2.49E-01
 
 ETA2
+       -3.80E-02  1.63E-01
 
 ETA3
+        4.71E-02 -1.52E-02  8.69E-02
 
 ETA4
+        3.18E-02  5.11E-02 -2.08E-02  2.24E-01
 
 ETA5
+        2.73E-02  2.64E-02 -2.32E-03 -3.34E-02  1.69E-01
 
 ETA6
+       -2.74E-02  8.22E-03  2.64E-02  1.84E-02 -8.02E-02  1.90E-01
 
 ETA7
+        2.97E-02 -3.96E-02  3.17E-02 -7.25E-02  2.44E-02  4.85E-03  2.13E-01
 
 ETA8
+        9.81E-02  8.20E-02  3.57E-02  4.41E-02  9.85E-04 -5.14E-02  5.57E-02  1.99E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1      EPS2     
 
 EPS1
+        9.27E-03
 
 EPS2
+        0.00E+00  2.24E-02
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        4.99E-01
 
 ETA2
+       -1.89E-01  4.03E-01
 
 ETA3
+        3.21E-01 -1.28E-01  2.95E-01
 
 ETA4
+        1.35E-01  2.68E-01 -1.49E-01  4.73E-01
 
 ETA5
+        1.33E-01  1.59E-01 -1.91E-02 -1.72E-01  4.12E-01
 
 ETA6
+       -1.26E-01  4.67E-02  2.05E-01  8.92E-02 -4.47E-01  4.36E-01
 
 ETA7
+        1.29E-01 -2.12E-01  2.33E-01 -3.32E-01  1.28E-01  2.41E-02  4.62E-01
 
 ETA8
+        4.42E-01  4.56E-01  2.72E-01  2.09E-01  5.37E-03 -2.64E-01  2.70E-01  4.46E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1      EPS2     
 
 EPS1
+        9.63E-02
 
 EPS2
+        0.00E+00  1.50E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                          STANDARD ERROR OF ESTIMATE (S)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8     
 
         3.19E-01  2.38E-01  2.27E-01  2.96E-01  2.05E-01  2.63E-01  1.50E-01  3.74E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        1.51E-01
 
 ETA2
+        2.63E-01  3.85E-01
 
 ETA3
+        1.10E-01  1.56E-01  1.93E-01
 
 ETA4
+        9.06E-02  1.54E-01  1.00E-01  1.71E-01
 
 ETA5
+        1.48E-01  9.39E-02  8.81E-02  1.58E-01  3.04E-01
 
 ETA6
+        1.14E-01  1.57E-01  1.38E-01  1.62E-01  7.80E-02  1.51E-01
 
 ETA7
+        1.83E-01  1.24E-01  1.15E-01  9.59E-02  1.17E-01  1.04E-01  1.51E-01
 
 ETA8
+        1.99E-01  1.40E-01  6.45E-02  9.48E-02  1.57E-01  1.94E-01  1.32E-01  2.40E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1      EPS2     
 
 EPS1
+        3.23E-03
 
 EPS2
+        0.00E+00  4.48E-03
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        1.51E-01
 
 ETA2
+        1.12E+00  4.77E-01
 
 ETA3
+        7.03E-01  1.12E+00  3.27E-01
 
 ETA4
+        3.75E-01  6.59E-01  8.12E-01  1.81E-01
 
 ETA5
+        7.47E-01  5.93E-01  7.27E-01  8.95E-01  3.69E-01
 
 ETA6
+        5.38E-01  8.88E-01  1.05E+00  7.73E-01  7.33E-01  1.73E-01
 
 ETA7
+        7.83E-01  5.37E-01  9.07E-01  4.14E-01  6.98E-01  5.18E-01  1.64E-01
 
 ETA8
+        6.23E-01  8.71E-01  3.37E-01  3.85E-01  8.59E-01  8.93E-01  5.82E-01  2.69E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1      EPS2     
 
 EPS1
+        1.68E-02
 
 EPS2
+       .........  1.49E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (S)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 TH 1
+        1.02E-01
 
 TH 2
+       -1.90E-02  5.67E-02
 
 TH 3
+        1.90E-02 -7.01E-04  5.16E-02
 
 TH 4
+        6.47E-02 -2.23E-02 -1.67E-02  8.77E-02
 
 TH 5
+       -4.43E-02  2.38E-02  5.84E-03 -4.36E-02  4.19E-02
 
 TH 6
+       -4.78E-02  6.00E-03 -1.82E-02 -2.23E-02  2.31E-02  6.90E-02
 
 TH 7
+       -1.33E-02 -5.87E-03  9.56E-03 -2.02E-02  1.07E-02  9.99E-03  2.25E-02
 
 TH 8
+        8.67E-02 -3.79E-02 -1.06E-03  8.62E-02 -5.28E-02 -3.56E-02 -9.09E-03  1.40E-01
 
 OM11
+       -2.32E-02 -2.20E-03 -2.09E-02 -6.74E-04  7.60E-03  1.51E-02 -1.75E-04 -6.91E-03  2.28E-02
 
 OM12
+       -1.51E-02  4.13E-02  1.72E-02 -3.49E-02  2.47E-02  3.20E-03 -2.66E-03 -6.43E-02 -1.28E-02  6.90E-02
 
 OM13
+       -1.17E-02 -1.24E-02 -7.77E-04 -4.96E-03  2.37E-03  6.38E-03  4.11E-03  2.70E-03  8.52E-03 -1.78E-02  1.21E-02
 
 OM14
+       -8.05E-03 -1.89E-03  5.35E-03 -1.34E-02  4.58E-03 -5.24E-03  6.73E-03 -1.12E-02 -7.29E-04 -5.05E-04  2.29E-03  8.20E-03
 
 OM15
+        1.04E-02 -1.01E-02 -8.40E-03  1.84E-02 -7.55E-03  1.29E-02 -6.03E-03  1.37E-02  3.38E-03 -5.62E-03  6.40E-04 -9.52E-03
          2.20E-02
 
 OM16
+       -1.26E-02  1.60E-02  3.07E-03 -1.67E-02  1.19E-02 -1.46E-06  1.51E-03 -1.42E-02  4.83E-04  1.23E-02 -1.49E-03  2.55E-03
         -1.04E-02  1.31E-02
 
 OM17
+        5.76E-03 -1.13E-02 -2.71E-02  3.15E-02 -1.08E-02  6.15E-03 -7.71E-03  3.64E-02  1.78E-02 -3.10E-02  6.29E-03 -7.04E-03
          1.08E-02 -5.82E-03  3.35E-02
 
 OM18
+       -1.96E-02  5.85E-03 -3.45E-02  1.14E-02  5.52E-03  2.65E-02 -6.20E-03 -1.45E-03  2.21E-02 -5.82E-03  1.04E-03 -7.63E-03
          1.11E-02 -1.99E-03  2.68E-02  3.97E-02
 
 OM22
+       -2.34E-02 -5.30E-02 -5.31E-02  2.63E-02 -2.16E-02  3.37E-02  1.68E-03  4.65E-02  3.07E-02 -8.30E-02  2.54E-02 -3.43E-03
          1.83E-02 -1.77E-02  4.54E-02  3.15E-02  1.48E-01
 
 OM23
+        6.72E-03  1.47E-02  2.21E-02 -1.75E-02  8.42E-03 -1.32E-02  4.48E-03 -2.16E-02 -1.40E-02  2.77E-02 -9.94E-03  5.02E-03
         -1.08E-02  5.87E-03 -2.13E-02 -1.71E-02 -5.05E-02  2.44E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM24
+       -1.06E-02 -1.25E-02 -1.76E-02  9.05E-03 -5.20E-03  1.76E-02 -5.67E-03  3.66E-03  8.64E-03 -8.74E-03  3.38E-03 -7.90E-03
          1.62E-02 -7.10E-03  1.23E-02  1.63E-02  3.57E-02 -1.64E-02  2.39E-02
 
 OM25
+        5.36E-03  1.05E-03  1.04E-03 -1.46E-03 -1.21E-03 -4.71E-03  3.11E-03  2.35E-03 -9.49E-04 -2.69E-03 -4.35E-04  3.01E-03
         -5.46E-03  2.20E-03 -1.84E-03 -5.72E-03 -2.69E-03  4.31E-03 -8.45E-03  8.82E-03
 
 OM26
+       -1.55E-02  2.01E-02 -1.52E-02 -7.15E-03  1.28E-02  2.27E-02 -3.48E-03 -1.94E-02  7.76E-03  1.73E-02 -5.53E-03 -6.70E-03
          5.72E-03  5.87E-03  3.95E-03  1.86E-02 -5.62E-03 -8.25E-04  7.35E-03 -2.52E-03  2.45E-02
 
 OM27
+       -8.86E-03  2.17E-02  4.24E-03 -1.60E-02  1.42E-02 -1.64E-03  1.08E-03 -2.63E-02 -7.50E-04  2.35E-02 -6.11E-03  1.90E-03
         -8.15E-03  8.82E-03 -9.24E-03 -4.11E-04 -3.44E-02  1.24E-02 -1.02E-02  2.90E-03  8.59E-03  1.54E-02
 
 OM28
+       -3.08E-02  1.32E-02 -9.92E-03 -2.34E-02  1.50E-02  1.31E-02  6.59E-04 -4.29E-02  7.50E-03  1.99E-02 -1.14E-03  9.80E-04
         -1.45E-03  5.52E-03 -6.27E-03  7.20E-03 -5.30E-03  1.59E-03  5.01E-03 -7.19E-04  1.04E-02  8.33E-03  1.95E-02
 
 OM33
+        3.75E-02 -2.89E-02 -8.97E-03  4.01E-02 -3.03E-02 -1.35E-02 -6.20E-03  6.35E-02 -2.24E-04 -4.05E-02  5.03E-03 -5.06E-03
          1.08E-02 -1.17E-02  2.07E-02  1.44E-03  4.43E-02 -1.62E-02  8.17E-03  4.52E-04 -1.11E-02 -1.81E-02 -1.92E-02  3.73E-02
 
 OM34
+        8.68E-03 -1.23E-02  5.74E-03  7.57E-03 -6.42E-03 -1.25E-02  4.33E-03  1.81E-02 -1.13E-03 -1.52E-02  4.17E-03  4.09E-03
         -4.81E-03 -1.71E-03  2.32E-03 -7.97E-03  7.50E-03 -4.72E-04 -5.59E-03  3.55E-03 -1.15E-02 -4.04E-03 -7.48E-03  8.24E-03
         1.00E-02
 
 OM35
+        2.22E-03  3.07E-03 -9.26E-03  6.24E-03 -2.22E-03  9.86E-03 -4.55E-03 -1.78E-04  1.58E-03  5.11E-04 -2.28E-03 -4.04E-03
          7.84E-03 -3.61E-03  4.79E-03  9.09E-03  7.31E-03 -3.75E-03  5.65E-03 -2.46E-03  6.41E-03 -1.68E-03  5.47E-04  2.50E-03
        -5.83E-03  7.77E-03
 
 OM36
+        1.69E-02  3.63E-03 -1.13E-02  1.61E-02 -1.31E-02 -1.71E-02 -1.18E-02  1.56E-02 -2.59E-04  1.26E-03 -6.53E-03 -4.64E-03
          2.44E-03  1.82E-03  7.05E-03  4.10E-03 -2.42E-03 -1.13E-03  3.74E-03  2.34E-04  4.81E-03  1.37E-03  6.75E-04  8.04E-03
        -1.69E-03  1.05E-03  1.89E-02
 
 OM37
+       -1.62E-02  7.68E-03 -1.61E-02 -3.61E-03  4.02E-03  1.34E-02 -5.59E-03 -1.45E-02  9.55E-03  4.41E-03  9.96E-04 -4.23E-03
          6.58E-03  1.29E-03  6.10E-03  1.42E-02  1.16E-02 -8.79E-03  1.07E-02 -4.39E-03  1.13E-02  8.83E-04  9.80E-03 -3.89E-03
        -8.31E-03  5.50E-03  4.17E-03  1.33E-02
 
 OM38
+        8.01E-03 -9.41E-03 -2.03E-03  9.03E-03 -6.87E-03  7.26E-05 -8.78E-04  1.51E-02  1.56E-03 -1.24E-02  3.46E-03 -1.42E-03
          4.41E-03 -3.55E-03  6.19E-03  1.30E-03  1.45E-02 -5.62E-03  2.85E-03 -2.50E-04 -2.38E-03 -5.12E-03 -4.94E-03  1.01E-02
         1.64E-03  1.52E-03  3.13E-04  2.99E-04  4.17E-03
 
 OM44
+       -2.55E-02 -1.15E-02  1.42E-02 -2.49E-02  1.12E-02  3.24E-03  8.49E-03 -2.13E-02  2.69E-03 -1.49E-03  1.07E-02  5.33E-03
         -1.79E-03 -8.21E-04 -8.04E-03 -8.77E-03  6.22E-03 -1.31E-03  4.82E-03 -3.94E-03 -9.69E-03 -3.45E-03  5.53E-03 -7.77E-03
         3.27E-03 -6.89E-03 -9.94E-03 -8.43E-04 -1.02E-03  2.92E-02
 
 OM45
+       -1.41E-02  2.29E-02 -6.04E-03 -1.35E-02  1.41E-02  1.95E-02 -2.15E-03 -3.38E-02  1.21E-03  2.92E-02 -7.38E-03 -3.21E-03
          3.32E-03  5.02E-03 -7.82E-03  9.40E-03 -2.42E-02  5.60E-03  1.26E-03 -1.58E-03  1.91E-02  1.09E-02  1.21E-02 -1.88E-02
        -1.27E-02  6.50E-03 -2.84E-04  9.72E-03 -3.84E-03 -8.47E-03  2.49E-02
 
 OM46
+       -1.95E-02 -1.17E-02 -4.87E-03 -9.87E-03  6.13E-03  1.67E-02  1.23E-02  4.94E-03  9.47E-03 -2.49E-02  1.21E-02  4.31E-03
         -5.05E-03  2.33E-03  5.83E-03  7.00E-04  3.72E-02 -9.95E-03 -1.72E-04  5.06E-03 -2.80E-03 -5.41E-03 -1.36E-03  5.03E-03
         6.22E-03 -3.58E-03 -1.07E-02 -1.74E-03  3.05E-03  6.67E-03 -7.97E-03  2.63E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM47
+       -9.81E-03  8.20E-03 -1.82E-03 -5.26E-03  8.93E-03  9.17E-03  5.90E-03 -6.45E-03  4.75E-03  1.25E-03  1.14E-03  1.92E-03
         -4.90E-03  4.51E-03  9.08E-04  4.55E-03 -2.87E-03  1.25E-04 -4.99E-03  1.85E-03  4.17E-03  4.97E-03  2.29E-03 -7.02E-03
        -1.02E-04 -5.38E-04 -5.10E-03  6.30E-04 -1.24E-03 -3.09E-03  4.72E-03  6.54E-03  9.20E-03
 
 OM48
+       -1.80E-02 -3.07E-03 -1.30E-02 -3.10E-03  3.31E-03  1.07E-02  1.49E-03 -1.01E-02  8.36E-03 -6.03E-03  3.70E-03  6.90E-04
          2.71E-03 -2.62E-03  6.12E-03  1.04E-02  2.10E-02 -8.15E-03  8.41E-03 -3.03E-03  2.44E-03 -2.24E-03  5.87E-03 -1.37E-03
        -1.26E-03  1.82E-03 -1.66E-03  5.39E-03 -4.69E-05  4.96E-03  5.70E-04  3.39E-03  1.34E-03  8.99E-03
 
 OM55
+       -4.79E-02  3.57E-02  3.23E-03 -5.48E-02  3.66E-02  1.56E-02  6.78E-03 -9.24E-02  6.83E-03  5.94E-02 -6.96E-03  3.48E-03
         -7.39E-03  1.51E-02 -2.60E-02  1.48E-03 -5.63E-02  1.96E-02 -3.43E-03  2.49E-03  2.10E-02  2.69E-02  3.61E-02 -4.93E-02
        -1.35E-02 -4.16E-03 -4.34E-04  1.18E-02 -1.34E-02  1.23E-02  2.79E-02 -9.63E-03  6.04E-03  3.19E-03  9.24E-02
 
 OM56
+       -3.61E-03  9.63E-03 -8.32E-03 -2.27E-03  1.45E-03  2.56E-03 -5.43E-03 -8.71E-03  1.57E-03  7.46E-03 -3.65E-03 -2.12E-03
          3.33E-04  2.85E-03  1.34E-03  6.12E-03 -3.41E-03  8.30E-05  1.85E-03 -5.45E-06  7.33E-03  4.06E-03  4.87E-03 -3.71E-03
        -4.99E-03  2.91E-03  5.62E-03  5.09E-03 -1.17E-03 -5.48E-03  6.60E-03 -4.58E-03 -1.95E-04  7.92E-04  7.78E-03  6.08E-03
 
 OM57
+        1.96E-03 -7.52E-03 -1.84E-02  1.63E-02 -5.70E-03  4.91E-03 -3.86E-03  1.99E-02  1.02E-02 -2.07E-02  4.77E-03 -3.29E-03
          5.61E-03 -3.41E-03  1.79E-02  1.41E-02  3.20E-02 -1.44E-02  6.61E-03  2.04E-04  2.05E-03 -5.27E-03 -3.32E-03  1.32E-02
         1.16E-03  3.34E-03  2.90E-03  4.35E-03  4.49E-03 -5.44E-03 -3.85E-03  6.83E-03  7.83E-04  4.17E-03 -1.63E-02  1.36E-03
          1.37E-02
 
 OM58
+        2.29E-02 -1.76E-02 -1.57E-02  3.35E-02 -1.89E-02 -1.96E-03 -5.01E-03  4.41E-02  6.44E-03 -3.15E-02  4.52E-03 -4.99E-03
          1.14E-02 -9.09E-03  2.19E-02  9.40E-03  3.87E-02 -1.56E-02  6.78E-03  3.30E-03 -3.99E-03 -1.17E-02 -1.10E-02  2.56E-02
         5.37E-03  3.79E-03  5.23E-03 -9.02E-05  7.46E-03 -9.07E-03 -1.01E-02  6.32E-03 -2.01E-03  1.14E-03 -3.15E-02 -1.74E-03
          1.40E-02  2.46E-02
 
 OM66
+       -2.69E-02  9.32E-03 -6.98E-03 -2.07E-02  1.44E-02  2.51E-02  4.20E-03 -2.67E-02  3.38E-03  8.50E-03  5.48E-04 -6.10E-05
          2.27E-03  2.20E-03 -4.50E-03  8.25E-03  7.87E-03 -1.29E-03  6.36E-03 -5.23E-03  1.16E-02  1.10E-03  7.94E-03 -9.98E-03
        -8.83E-03  6.08E-03 -6.38E-03  8.43E-03 -9.18E-04  1.41E-03  1.24E-02  4.43E-03  3.05E-03  4.73E-03  7.91E-03  1.23E-03
         -1.44E-03 -7.45E-03  2.28E-02
 
 OM67
+        1.09E-02 -1.77E-02  2.70E-03  8.46E-03 -1.28E-02 -7.52E-03 -3.57E-04  1.56E-02 -4.27E-03 -1.18E-02  2.19E-03  1.51E-04
          3.70E-03 -6.40E-03  8.18E-05 -7.16E-03  1.49E-02 -2.86E-03  4.65E-03 -1.03E-03 -8.85E-03 -9.15E-03 -5.48E-03  1.17E-02
         4.02E-03 -1.05E-03  1.83E-03 -3.35E-03  2.59E-03  4.36E-03 -9.22E-03 -1.18E-04 -6.56E-03 -2.45E-04 -1.44E-02 -3.36E-03
         -5.34E-04  4.91E-03 -2.92E-03  1.08E-02
 
 OM68
+        5.13E-02 -8.00E-03  3.66E-03  4.04E-02 -2.67E-02 -2.62E-02 -1.15E-02  5.72E-02 -9.87E-03 -1.07E-02 -7.46E-03 -8.05E-03
          6.81E-03 -3.09E-03  1.00E-02 -4.51E-03 -5.78E-03 -2.44E-04  1.25E-04  5.27E-04 -2.49E-03 -6.25E-03 -1.69E-02  2.44E-02
         4.30E-03  6.95E-04  1.67E-02 -6.14E-03  4.29E-03 -1.57E-02 -9.42E-03 -1.08E-02 -7.14E-03 -9.66E-03 -3.02E-02 -3.48E-04
          3.02E-03  1.43E-02 -1.49E-02  7.10E-03  3.75E-02
 
 OM77
+       -5.14E-03  1.33E-02  2.07E-03 -5.35E-03  5.72E-03 -4.70E-03 -7.74E-03 -2.65E-02 -2.73E-03  2.54E-02 -7.93E-03 -1.31E-03
          1.30E-03 -6.77E-04 -7.04E-03  2.04E-03 -2.87E-02  9.31E-03  8.88E-04 -2.86E-03  3.57E-03  6.49E-03  8.62E-03 -1.40E-02
        -5.50E-03  1.39E-03  3.51E-03  2.94E-03 -5.36E-03  2.12E-03  7.77E-03 -1.69E-02 -4.63E-03  5.26E-04  2.21E-02  3.57E-03
         -6.96E-03 -9.62E-03  6.82E-05 -2.08E-03 -3.90E-03  2.29E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM78
+        1.50E-02  4.60E-03 -1.05E-02  2.45E-02 -5.65E-03  1.88E-03 -9.37E-03  1.86E-02  4.63E-03 -3.73E-04 -4.86E-03 -7.56E-03
          8.73E-03 -3.00E-03  1.52E-02  1.70E-02 -8.93E-04 -4.55E-03  4.85E-03 -2.67E-03  7.88E-03  4.26E-04 -3.82E-03  6.54E-03
        -3.50E-03  5.82E-03  6.37E-03  3.70E-03  1.82E-03 -1.24E-02  4.02E-03 -8.32E-03  9.50E-04 -3.02E-04 -9.04E-03  3.03E-03
          6.50E-03  8.93E-03 -2.50E-03 -4.12E-03  1.10E-02  4.54E-03  1.75E-02
 
 OM88
+       -4.07E-02  2.82E-02 -2.76E-02 -1.54E-02  2.60E-02  3.67E-02 -3.42E-03 -4.78E-02  1.98E-02  2.64E-02 -4.63E-03 -5.61E-03
          6.39E-03  5.84E-03  9.77E-03  3.68E-02 -4.83E-03 -3.81E-03  1.05E-02 -5.26E-03  2.92E-02  1.41E-02  2.19E-02 -2.47E-02
        -1.66E-02  9.81E-03 -6.46E-05  1.90E-02 -4.55E-03 -6.82E-03  2.74E-02 -5.25E-03  9.32E-03  1.02E-02  3.81E-02  1.03E-02
          4.89E-03 -7.82E-03  1.74E-02 -1.66E-02 -2.09E-02  1.25E-02  1.41E-02  5.77E-02
 
 SG11
+        7.77E-04 -4.15E-04 -4.96E-05  8.12E-04 -5.10E-04 -3.14E-04 -1.56E-04  1.03E-03 -6.71E-05 -4.59E-04  2.84E-06 -1.19E-04
          1.96E-04 -1.78E-04  2.70E-04 -1.99E-05  3.52E-04 -1.99E-04  9.46E-05  7.30E-06 -1.75E-04 -2.25E-04 -3.00E-04  5.18E-04
         1.21E-04  2.43E-05  1.52E-04 -8.44E-05  1.31E-04 -1.78E-04 -2.34E-04 -3.10E-05 -9.46E-05 -7.05E-05 -6.52E-04 -5.92E-05
          1.71E-04  3.73E-04 -2.64E-04  1.46E-04  4.59E-04 -1.61E-04  1.65E-04 -3.65E-04  1.04E-05
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG22
+       -1.85E-06 -7.04E-04  1.97E-05  9.38E-05 -2.24E-04  3.03E-05  1.86E-04  6.61E-04  6.08E-05 -9.32E-04  3.41E-04  9.31E-05
          2.57E-05 -2.03E-04  2.27E-04 -1.06E-04  1.14E-03 -3.39E-04  1.09E-04 -2.01E-05 -4.02E-04 -3.86E-04 -3.08E-04  4.99E-04
         2.37E-04 -5.48E-05 -2.29E-04 -1.30E-04  1.78E-04  2.64E-04 -4.81E-04  4.22E-04 -4.50E-05  8.91E-05 -8.79E-04 -1.77E-04
          1.99E-04  2.94E-04 -1.42E-05  2.14E-04 -6.08E-05 -3.71E-04 -1.58E-04 -5.34E-04  4.06E-06  0.00E+00  2.00E-05
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (S)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 TH 1
+        3.19E-01
 
 TH 2
+       -2.49E-01  2.38E-01
 
 TH 3
+        2.61E-01 -1.30E-02  2.27E-01
 
 TH 4
+        6.84E-01 -3.17E-01 -2.47E-01  2.96E-01
 
 TH 5
+       -6.77E-01  4.89E-01  1.26E-01 -7.19E-01  2.05E-01
 
 TH 6
+       -5.70E-01  9.59E-02 -3.05E-01 -2.87E-01  4.30E-01  2.63E-01
 
 TH 7
+       -2.78E-01 -1.64E-01  2.81E-01 -4.54E-01  3.50E-01  2.53E-01  1.50E-01
 
 TH 8
+        7.26E-01 -4.26E-01 -1.24E-02  7.79E-01 -6.91E-01 -3.63E-01 -1.62E-01  3.74E-01
 
 OM11
+       -4.81E-01 -6.11E-02 -6.10E-01 -1.51E-02  2.46E-01  3.82E-01 -7.73E-03 -1.22E-01  1.51E-01
 
 OM12
+       -1.80E-01  6.61E-01  2.89E-01 -4.49E-01  4.60E-01  4.63E-02 -6.74E-02 -6.55E-01 -3.22E-01  2.63E-01
 
 OM13
+       -3.32E-01 -4.74E-01 -3.11E-02 -1.52E-01  1.05E-01  2.21E-01  2.50E-01  6.59E-02  5.13E-01 -6.18E-01  1.10E-01
 
 OM14
+       -2.78E-01 -8.75E-02  2.60E-01 -5.00E-01  2.47E-01 -2.20E-01  4.95E-01 -3.30E-01 -5.33E-02 -2.12E-02  2.30E-01  9.06E-02
 
 OM15
+        2.19E-01 -2.86E-01 -2.49E-01  4.20E-01 -2.49E-01  3.31E-01 -2.71E-01  2.48E-01  1.51E-01 -1.44E-01  3.93E-02 -7.09E-01
          1.48E-01
 
 OM16
+       -3.46E-01  5.90E-01  1.18E-01 -4.93E-01  5.10E-01 -4.85E-05  8.79E-02 -3.32E-01  2.80E-02  4.09E-01 -1.19E-01  2.47E-01
         -6.14E-01  1.14E-01
 
 OM17
+        9.85E-02 -2.60E-01 -6.52E-01  5.81E-01 -2.88E-01  1.28E-01 -2.81E-01  5.33E-01  6.45E-01 -6.44E-01  3.13E-01 -4.25E-01
          3.99E-01 -2.78E-01  1.83E-01
 
 OM18
+       -3.08E-01  1.23E-01 -7.61E-01  1.94E-01  1.35E-01  5.06E-01 -2.07E-01 -1.95E-02  7.33E-01 -1.11E-01  4.73E-02 -4.22E-01
          3.74E-01 -8.72E-02  7.34E-01  1.99E-01
 
 OM22
+       -1.91E-01 -5.79E-01 -6.07E-01  2.30E-01 -2.74E-01  3.33E-01  2.92E-02  3.24E-01  5.28E-01 -8.21E-01  6.02E-01 -9.84E-02
          3.20E-01 -4.02E-01  6.45E-01  4.10E-01  3.85E-01
 
 OM23
+        1.35E-01  3.97E-01  6.24E-01 -3.80E-01  2.63E-01 -3.22E-01  1.91E-01 -3.70E-01 -5.95E-01  6.75E-01 -5.80E-01  3.55E-01
         -4.67E-01  3.29E-01 -7.46E-01 -5.50E-01 -8.42E-01  1.56E-01
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM24
+       -2.14E-01 -3.41E-01 -5.03E-01  1.98E-01 -1.64E-01  4.34E-01 -2.45E-01  6.35E-02  3.70E-01 -2.16E-01  1.99E-01 -5.64E-01
          7.06E-01 -4.02E-01  4.36E-01  5.30E-01  6.00E-01 -6.81E-01  1.54E-01
 
 OM25
+        1.79E-01  4.68E-02  4.87E-02 -5.26E-02 -6.28E-02 -1.91E-01  2.20E-01  6.70E-02 -6.69E-02 -1.09E-01 -4.21E-02  3.53E-01
         -3.92E-01  2.05E-01 -1.07E-01 -3.05E-01 -7.44E-02  2.94E-01 -5.82E-01  9.39E-02
 
 OM26
+       -3.09E-01  5.39E-01 -4.28E-01 -1.54E-01  3.98E-01  5.52E-01 -1.48E-01 -3.31E-01  3.28E-01  4.20E-01 -3.22E-01 -4.72E-01
          2.46E-01  3.28E-01  1.38E-01  5.97E-01 -9.33E-02 -3.38E-02  3.04E-01 -1.71E-01  1.57E-01
 
 OM27
+       -2.23E-01  7.34E-01  1.50E-01 -4.36E-01  5.57E-01 -5.02E-02  5.81E-02 -5.66E-01 -4.00E-02  7.19E-01 -4.47E-01  1.68E-01
         -4.43E-01  6.21E-01 -4.06E-01 -1.66E-02 -7.19E-01  6.39E-01 -5.32E-01  2.48E-01  4.41E-01  1.24E-01
 
 OM28
+       -6.90E-01  3.98E-01 -3.12E-01 -5.65E-01  5.26E-01  3.56E-01  3.14E-02 -8.22E-01  3.55E-01  5.43E-01 -7.43E-02  7.74E-02
         -7.01E-02  3.45E-01 -2.45E-01  2.58E-01 -9.86E-02  7.31E-02  2.32E-01 -5.48E-02  4.78E-01  4.80E-01  1.40E-01
 
 OM33
+        6.09E-01 -6.30E-01 -2.04E-01  7.02E-01 -7.67E-01 -2.67E-01 -2.14E-01  8.80E-01 -7.67E-03 -7.99E-01  2.37E-01 -2.89E-01
          3.77E-01 -5.31E-01  5.86E-01  3.75E-02  5.97E-01 -5.38E-01  2.74E-01  2.49E-02 -3.66E-01 -7.53E-01 -7.11E-01  1.93E-01
 
 OM34
+        2.71E-01 -5.17E-01  2.52E-01  2.55E-01 -3.13E-01 -4.76E-01  2.88E-01  4.83E-01 -7.49E-02 -5.78E-01  3.79E-01  4.51E-01
         -3.24E-01 -1.49E-01  1.26E-01 -3.99E-01  1.95E-01 -3.02E-02 -3.61E-01  3.77E-01 -7.34E-01 -3.25E-01 -5.34E-01  4.26E-01
         1.00E-01
 
 OM35
+        7.90E-02  1.46E-01 -4.62E-01  2.39E-01 -1.23E-01  4.26E-01 -3.44E-01 -5.42E-03  1.19E-01  2.21E-02 -2.36E-01 -5.06E-01
          6.00E-01 -3.58E-01  2.97E-01  5.17E-01  2.15E-01 -2.73E-01  4.15E-01 -2.98E-01  4.65E-01 -1.53E-01  4.44E-02  1.47E-01
        -6.61E-01  8.81E-02
 
 OM36
+        3.86E-01  1.11E-01 -3.61E-01  3.96E-01 -4.67E-01 -4.73E-01 -5.70E-01  3.03E-01 -1.24E-02  3.48E-02 -4.32E-01 -3.72E-01
          1.20E-01  1.16E-01  2.80E-01  1.49E-01 -4.58E-02 -5.27E-02  1.76E-01  1.81E-02  2.23E-01  8.01E-02  3.51E-02  3.03E-01
        -1.23E-01  8.63E-02  1.38E-01
 
 OM37
+       -4.40E-01  2.80E-01 -6.13E-01 -1.06E-01  1.71E-01  4.44E-01 -3.24E-01 -3.38E-01  5.49E-01  1.46E-01  7.86E-02 -4.05E-01
          3.85E-01  9.78E-02  2.89E-01  6.20E-01  2.62E-01 -4.89E-01  6.01E-01 -4.06E-01  6.25E-01  6.17E-02  6.09E-01 -1.75E-01
        -7.20E-01  5.41E-01  2.63E-01  1.15E-01
 
 OM38
+        3.89E-01 -6.13E-01 -1.38E-01  4.73E-01 -5.20E-01  4.28E-03 -9.07E-02  6.27E-01  1.60E-01 -7.29E-01  4.88E-01 -2.42E-01
          4.60E-01 -4.81E-01  5.24E-01  1.01E-01  5.84E-01 -5.58E-01  2.86E-01 -4.13E-02 -2.35E-01 -6.39E-01 -5.48E-01  8.12E-01
         2.53E-01  2.67E-01  3.53E-02  4.02E-02  6.45E-02
 
 OM44
+       -4.68E-01 -2.82E-01  3.66E-01 -4.92E-01  3.20E-01  7.21E-02  3.31E-01 -3.34E-01  1.04E-01 -3.31E-02  5.69E-01  3.45E-01
         -7.05E-02 -4.21E-02 -2.57E-01 -2.58E-01  9.47E-02 -4.92E-02  1.83E-01 -2.46E-01 -3.62E-01 -1.62E-01  2.32E-01 -2.36E-01
         1.91E-01 -4.58E-01 -4.23E-01 -4.28E-02 -9.24E-02  1.71E-01
 
 OM45
+       -2.81E-01  6.09E-01 -1.69E-01 -2.90E-01  4.36E-01  4.72E-01 -9.08E-02 -5.74E-01  5.09E-02  7.04E-01 -4.26E-01 -2.25E-01
          1.42E-01  2.78E-01 -2.71E-01  2.99E-01 -3.99E-01  2.27E-01  5.16E-02 -1.07E-01  7.72E-01  5.57E-01  5.50E-01 -6.18E-01
        -8.05E-01  4.67E-01 -1.31E-02  5.35E-01 -3.78E-01 -3.14E-01  1.58E-01
 
 OM46
+       -3.77E-01 -3.04E-01 -1.32E-01 -2.06E-01  1.85E-01  3.91E-01  5.05E-01  8.16E-02  3.87E-01 -5.84E-01  6.76E-01  2.94E-01
         -2.10E-01  1.26E-01  1.97E-01  2.16E-02  5.96E-01 -3.93E-01 -6.86E-03  3.32E-01 -1.10E-01 -2.68E-01 -6.02E-02  1.61E-01
         3.83E-01 -2.51E-01 -4.79E-01 -9.31E-02  2.91E-01  2.41E-01 -3.12E-01  1.62E-01
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM47
+       -3.20E-01  3.59E-01 -8.34E-02 -1.85E-01  4.55E-01  3.64E-01  4.10E-01 -1.80E-01  3.28E-01  4.97E-02  1.08E-01  2.22E-01
         -3.44E-01  4.11E-01  5.17E-02  2.38E-01 -7.78E-02  8.35E-03 -3.37E-01  2.06E-01  2.78E-01  4.17E-01  1.71E-01 -3.79E-01
        -1.06E-02 -6.37E-02 -3.87E-01  5.70E-02 -2.00E-01 -1.89E-01  3.12E-01  4.21E-01  9.59E-02
 
 OM48
+       -5.96E-01 -1.36E-01 -6.03E-01 -1.10E-01  1.71E-01  4.30E-01  1.05E-01 -2.84E-01  5.84E-01 -2.42E-01  3.55E-01  8.04E-02
          1.93E-01 -2.42E-01  3.53E-01  5.49E-01  5.76E-01 -5.50E-01  5.74E-01 -3.41E-01  1.64E-01 -1.90E-01  4.43E-01 -7.50E-02
        -1.33E-01  2.18E-01 -1.27E-01  4.93E-01 -7.67E-03  3.06E-01  3.81E-02  2.20E-01  1.47E-01  9.48E-02
 
 OM55
+       -4.94E-01  4.94E-01  4.68E-02 -6.08E-01  5.88E-01  1.95E-01  1.49E-01 -8.13E-01  1.49E-01  7.45E-01 -2.08E-01  1.26E-01
         -1.64E-01  4.34E-01 -4.67E-01  2.44E-02 -4.82E-01  4.14E-01 -7.31E-02  8.73E-02  4.41E-01  7.12E-01  8.50E-01 -8.41E-01
        -4.44E-01 -1.55E-01 -1.04E-02  3.37E-01 -6.82E-01  2.37E-01  5.82E-01 -1.95E-01  2.07E-01  1.11E-01  3.04E-01
 
 OM56
+       -1.45E-01  5.19E-01 -4.70E-01 -9.84E-02  9.11E-02  1.25E-01 -4.64E-01 -2.99E-01  1.33E-01  3.64E-01 -4.26E-01 -3.01E-01
          2.88E-02  3.19E-01  9.37E-02  3.94E-01 -1.14E-01  6.82E-03  1.53E-01 -7.43E-04  6.00E-01  4.19E-01  4.46E-01 -2.46E-01
        -6.38E-01  4.23E-01  5.23E-01  5.66E-01 -2.33E-01 -4.11E-01  5.37E-01 -3.62E-01 -2.61E-02  1.07E-01  3.28E-01  7.80E-02
 
 OM57
+        5.23E-02 -2.70E-01 -6.90E-01  4.71E-01 -2.38E-01  1.60E-01 -2.20E-01  4.54E-01  5.76E-01 -6.73E-01  3.71E-01 -3.10E-01
          3.23E-01 -2.55E-01  8.36E-01  6.03E-01  7.10E-01 -7.89E-01  3.65E-01  1.86E-02  1.12E-01 -3.62E-01 -2.03E-01  5.86E-01
         9.90E-02  3.24E-01  1.80E-01  3.22E-01  5.94E-01 -2.72E-01 -2.09E-01  3.60E-01  6.98E-02  3.76E-01 -4.59E-01  1.49E-01
          1.17E-01
 
 OM58
+        4.57E-01 -4.73E-01 -4.42E-01  7.21E-01 -5.90E-01 -4.75E-02 -2.13E-01  7.53E-01  2.72E-01 -7.65E-01  2.62E-01 -3.51E-01
          4.91E-01 -5.07E-01  7.62E-01  3.01E-01  6.41E-01 -6.37E-01  2.80E-01  2.24E-01 -1.62E-01 -6.03E-01 -5.02E-01  8.46E-01
         3.42E-01  2.75E-01  2.43E-01 -4.99E-03  7.37E-01 -3.39E-01 -4.07E-01  2.49E-01 -1.34E-01  7.69E-02 -6.60E-01 -1.42E-01
          7.63E-01  1.57E-01
 
 OM66
+       -5.58E-01  2.59E-01 -2.03E-01 -4.62E-01  4.67E-01  6.33E-01  1.85E-01 -4.73E-01  1.48E-01  2.14E-01  3.30E-02 -4.46E-03
          1.02E-01  1.28E-01 -1.63E-01  2.74E-01  1.35E-01 -5.48E-02  2.73E-01 -3.69E-01  4.92E-01  5.84E-02  3.76E-01 -3.42E-01
        -5.83E-01  4.57E-01 -3.07E-01  4.84E-01 -9.41E-02  5.47E-02  5.20E-01  1.81E-01  2.11E-01  3.30E-01  1.72E-01  1.05E-01
         -8.12E-02 -3.14E-01  1.51E-01
 
 OM67
+        3.29E-01 -7.19E-01  1.15E-01  2.76E-01 -6.03E-01 -2.76E-01 -2.29E-02  4.03E-01 -2.72E-01 -4.35E-01  1.92E-01  1.61E-02
          2.41E-01 -5.40E-01  4.31E-03 -3.46E-01  3.74E-01 -1.76E-01  2.90E-01 -1.05E-01 -5.45E-01 -7.10E-01 -3.78E-01  5.86E-01
         3.87E-01 -1.15E-01  1.28E-01 -2.80E-01  3.86E-01  2.46E-01 -5.64E-01 -7.02E-03 -6.59E-01 -2.49E-02 -4.58E-01 -4.15E-01
         -4.40E-02  3.02E-01 -1.86E-01  1.04E-01
 
 OM68
+        8.29E-01 -1.74E-01  8.32E-02  7.03E-01 -6.74E-01 -5.15E-01 -3.96E-01  7.90E-01 -3.37E-01 -2.10E-01 -3.50E-01 -4.59E-01
          2.37E-01 -1.39E-01  2.83E-01 -1.17E-01 -7.75E-02 -8.06E-03  4.16E-03  2.89E-02 -8.21E-02 -2.60E-01 -6.25E-01  6.53E-01
         2.22E-01  4.07E-02  6.25E-01 -2.75E-01  3.43E-01 -4.75E-01 -3.08E-01 -3.44E-01 -3.84E-01 -5.26E-01 -5.13E-01 -2.31E-02
          1.33E-01  4.72E-01 -5.10E-01  3.54E-01  1.94E-01
 
 OM77
+       -1.06E-01  3.68E-01  6.03E-02 -1.19E-01  1.84E-01 -1.18E-01 -3.41E-01 -4.69E-01 -1.19E-01  6.40E-01 -4.77E-01 -9.56E-02
          5.81E-02 -3.91E-02 -2.54E-01  6.77E-02 -4.92E-01  3.94E-01  3.80E-02 -2.01E-01  1.50E-01  3.45E-01  4.07E-01 -4.80E-01
        -3.62E-01  1.04E-01  1.69E-01  1.69E-01 -5.48E-01  8.19E-02  3.25E-01 -6.90E-01 -3.19E-01  3.66E-02  4.80E-01  3.02E-01
         -3.93E-01 -4.05E-01  2.98E-03 -1.32E-01 -1.33E-01  1.51E-01
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM78
+        3.54E-01  1.46E-01 -3.49E-01  6.25E-01 -2.09E-01  5.42E-02 -4.73E-01  3.76E-01  2.32E-01 -1.08E-02 -3.35E-01 -6.32E-01
          4.46E-01 -1.99E-01  6.29E-01  6.46E-01 -1.76E-02 -2.20E-01  2.38E-01 -2.15E-01  3.81E-01  2.60E-02 -2.07E-01  2.56E-01
        -2.65E-01  5.00E-01  3.50E-01  2.43E-01  2.13E-01 -5.47E-01  1.93E-01 -3.88E-01  7.50E-02 -2.41E-02 -2.25E-01  2.94E-01
          4.20E-01  4.31E-01 -1.25E-01 -3.01E-01  4.31E-01  2.27E-01  1.32E-01
 
 OM88
+       -5.30E-01  4.93E-01 -5.07E-01 -2.16E-01  5.29E-01  5.82E-01 -9.49E-02 -5.32E-01  5.45E-01  4.19E-01 -1.76E-01 -2.58E-01
          1.80E-01  2.13E-01  2.22E-01  7.69E-01 -5.23E-02 -1.02E-01  2.84E-01 -2.33E-01  7.76E-01  4.71E-01  6.53E-01 -5.32E-01
        -6.90E-01  4.64E-01 -1.95E-03  6.88E-01 -2.94E-01 -1.66E-01  7.24E-01 -1.35E-01  4.05E-01  4.50E-01  5.23E-01  5.50E-01
          1.74E-01 -2.08E-01  4.79E-01 -6.67E-01 -4.48E-01  3.45E-01  4.43E-01  2.40E-01
 
 SG11
+        7.54E-01 -5.40E-01 -6.77E-02  8.49E-01 -7.72E-01 -3.71E-01 -3.23E-01  8.54E-01 -1.38E-01 -5.41E-01  8.02E-03 -4.08E-01
          4.09E-01 -4.82E-01  4.58E-01 -3.10E-02  2.84E-01 -3.94E-01  1.90E-01  2.41E-02 -3.47E-01 -5.60E-01 -6.65E-01  8.31E-01
         3.73E-01  8.55E-02  3.41E-01 -2.27E-01  6.30E-01 -3.22E-01 -4.60E-01 -5.91E-02 -3.06E-01 -2.30E-01 -6.65E-01 -2.35E-01
          4.53E-01  7.37E-01 -5.42E-01  4.37E-01  7.34E-01 -3.29E-01  3.87E-01 -4.71E-01  3.23E-03
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG22
+       -1.29E-03 -6.61E-01  1.94E-02  7.07E-02 -2.45E-01  2.57E-02  2.77E-01  3.95E-01  9.00E-02 -7.93E-01  6.93E-01  2.30E-01
          3.87E-02 -3.97E-01  2.78E-01 -1.19E-01  6.65E-01 -4.85E-01  1.57E-01 -4.77E-02 -5.74E-01 -6.93E-01 -4.92E-01  5.77E-01
         5.29E-01 -1.39E-01 -3.72E-01 -2.51E-01  6.17E-01  3.45E-01 -6.82E-01  5.82E-01 -1.05E-01  2.10E-01 -6.46E-01 -5.07E-01
          3.80E-01  4.19E-01 -2.10E-02  4.62E-01 -7.01E-02 -5.47E-01 -2.66E-01 -4.97E-01  2.81E-01  0.00E+00  4.48E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (S)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 TH 1
+        3.15E+02
 
 TH 2
+        1.70E+02  4.35E+02
 
 TH 3
+       -5.24E+01  2.05E+01  5.36E+02
 
 TH 4
+       -3.17E+01 -1.78E+01  5.84E+01  2.84E+02
 
 TH 5
+       -8.66E+01 -1.25E+02 -8.10E+00  2.34E+01  3.64E+02
 
 TH 6
+       -2.45E+01 -6.73E+01 -7.75E+01 -5.39E+01  1.18E+02  2.84E+02
 
 TH 7
+        6.15E+01  1.62E+02 -1.34E+01  1.11E+02 -9.18E+01 -8.62E+01  3.72E+02
 
 TH 8
+       -2.10E+02 -2.87E+02 -1.05E+02 -9.54E+01  1.28E+02  1.49E+02 -2.32E+02  5.69E+02
 
 OM11
+        8.90E+01  4.18E+01 -1.00E+02  2.20E+01  2.20E+01  1.15E+02  1.17E+02 -4.98E+01  1.02E+03
 
 OM12
+        1.16E+02  8.86E+01  4.29E+01  4.35E+02  8.00E+01  4.57E+01  3.74E+02 -3.62E+02  9.59E+02  3.64E+03
 
 OM13
+       -2.21E+02 -3.45E+01 -2.74E+02  5.58E+01 -1.51E+02  1.26E+01  1.84E+02  1.99E+02 -9.01E+02 -8.29E+02  3.52E+03
 
 OM14
+        5.59E+01  4.85E+02  1.08E+02  8.68E+01 -1.04E+02  8.18E+01  1.31E+02 -1.65E+02 -2.92E+02  4.02E+02  3.16E+02  2.68E+03
 
 OM15
+        6.84E+01  1.45E+02 -1.82E+02 -9.70E+01 -4.24E+02 -3.67E+02 -1.81E+01 -1.03E+02 -7.55E+02 -1.84E+03  1.08E+03 -2.10E+02
          2.86E+03
 
 OM16
+        2.28E+02  5.40E+01 -2.45E+01  9.12E+01 -3.37E+02 -3.17E+02  1.28E+02 -2.39E+02  1.65E+02 -3.81E+02 -5.83E+02 -8.00E+02
          8.71E+02  1.83E+03
 
 OM17
+        1.79E+02  2.08E+02  1.88E+02  1.15E+02  1.46E+01  1.70E+02  1.33E+02 -2.34E+02  5.91E+02  2.13E+03 -1.48E+03  1.12E+03
         -1.72E+03 -4.29E+02  3.03E+03
 
 OM18
+       -1.06E+02 -3.00E+02  1.38E+02 -1.22E+02 -4.71E+01 -2.40E+02 -3.43E+02  2.71E+02 -1.20E+03 -3.07E+03  1.01E+03 -9.63E+02
          2.12E+03  6.30E+02 -2.62E+03  4.33E+03
 
 OM22
+        1.14E+02  2.02E+02  1.91E+01  2.66E+02 -5.85E+01 -3.51E+01  3.03E+02 -3.46E+02  4.11E+02  2.20E+03 -2.20E+02  2.13E+02
         -8.48E+02 -1.93E+02  1.10E+03 -1.85E+03  1.94E+03
 
 OM23
+        1.36E+01 -1.24E+02 -1.67E+02 -9.44E+01 -1.63E+02  1.76E+02 -2.09E+02  2.51E+02 -5.63E+02 -1.03E+03  1.77E+03 -4.13E+02
          8.51E+02 -6.03E+01 -8.53E+02  1.51E+03 -3.16E+02  3.56E+03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM24
+        4.88E+02  6.53E+02 -1.01E+02 -5.40E+01 -9.28E+01  1.14E+02  2.03E+02 -4.98E+02  2.53E+02  4.75E+02 -4.24E+02  1.98E+03
         -8.63E+02 -3.08E+02  1.55E+03 -1.23E+03 -4.71E+01 -3.66E+02  4.42E+03
 
 OM25
+        8.63E+01  3.95E+01 -1.96E+02 -6.21E+01 -2.01E+02 -4.24E+02 -1.09E+01 -1.80E+02 -6.87E+02 -1.93E+03  6.06E+02 -8.59E+02
          2.94E+03  1.00E+03 -2.03E+03  2.43E+03 -1.16E+03 -1.61E+02 -7.18E+02  5.90E+03
 
 OM26
+        6.26E+00 -4.50E+01  1.13E+02  1.25E+02 -3.51E+02 -4.27E+02  2.21E+02 -4.76E+01 -4.88E+01 -4.14E+02 -2.67E+02 -3.08E+02
          9.66E+02  1.24E+03 -6.34E+02  7.60E+02 -4.48E+02 -6.09E+02 -1.10E+03  1.78E+03  3.30E+03
 
 OM27
+        2.60E+02  2.17E+02 -1.97E+02  1.76E+02  8.13E+01  3.23E+02  2.73E+02 -2.14E+02  6.90E+02  2.21E+03 -6.32E+02  1.53E+03
         -1.96E+03 -6.36E+02  2.33E+03 -2.59E+03  1.32E+03 -9.66E+02  2.81E+03 -2.56E+03 -1.46E+03  4.53E+03
 
 OM28
+       -4.07E+02 -5.09E+02  1.35E+02 -4.11E+02 -8.59E+01 -6.95E+01 -4.54E+02  9.43E+02 -1.04E+03 -4.47E+03  1.05E+03 -1.16E+03
          2.47E+03  8.33E+02 -2.86E+03  4.48E+03 -3.18E+03  1.39E+03 -2.32E+03  2.81E+03  2.21E+03 -3.83E+03  8.55E+03
 
 OM33
+       -1.18E+02 -6.16E+01  5.80E+01  1.55E+02 -7.59E+00  1.04E+01 -1.53E+00 -1.29E+02  8.55E+00 -1.26E+01  3.12E+02  1.01E+02
          1.12E+02  4.32E+02  1.55E+02 -4.66E+02 -2.92E+02 -8.07E+01  2.45E+02 -1.23E+02  1.30E+02  3.73E+02  1.05E+02  2.86E+03
 
 OM34
+        1.17E+02 -9.22E+01  3.80E+01 -2.58E+02  7.59E+01  1.05E+02 -1.87E+02  1.50E+02  3.03E+00 -4.78E+02 -1.46E+02 -7.73E+02
          1.73E+02  2.28E+02 -3.85E+02  1.84E+02 -3.42E+02  2.00E+02 -1.68E+02  9.17E+02  2.07E+02 -3.75E+02  8.49E+02  6.45E+02
         2.93E+03
 
 OM35
+       -1.77E+02 -1.27E+02  6.12E+01  9.65E+01  2.49E+02  6.19E+01  9.78E+01  2.98E+02  6.40E+02  8.97E+02 -9.25E+02  5.74E+01
         -1.22E+03  3.63E+02  7.98E+02 -1.32E+03 -3.72E+01 -1.89E+03  7.37E+02 -5.50E+02  6.73E+02  1.23E+03 -3.07E+02  8.32E+02
         9.82E+02  4.13E+03
 
 OM36
+       -2.67E+01  1.62E+02  9.43E+01  6.30E+01  1.06E+02  3.23E+01  1.58E+02 -1.06E+02  4.64E+00  1.50E+02  2.36E+02  9.26E+01
          2.17E+02 -1.34E+02 -3.25E+02 -3.12E+01  2.62E+01 -6.48E+02  8.81E+01  4.43E+02  8.04E+01 -2.42E+02 -1.43E+02 -8.73E+01
        -2.19E+02  1.09E+03  1.79E+03
 
 OM37
+        1.76E+02 -2.65E+02 -2.43E+01 -1.58E+02  1.18E+02  1.79E+02 -2.56E+02  7.39E+01 -7.28E+02 -7.10E+02  1.03E+03 -3.57E+02
          7.54E+02 -4.15E+02 -6.83E+02  1.25E+03 -3.12E+02  2.48E+03 -2.77E+02  1.14E+03 -5.11E+02 -5.96E+02  5.43E+02 -6.13E+01
         1.34E+03 -1.84E+03 -9.20E+02  4.25E+03
 
 OM38
+        1.56E+02  2.60E+02 -1.67E+02  1.20E+02  3.25E+02 -1.28E+02  5.29E+01 -1.39E+02  9.33E+02  1.49E+03 -3.63E+03 -1.21E+02
         -1.36E+03  3.35E+02  1.37E+03 -8.28E+02  7.66E+02 -2.97E+03  6.51E+02 -2.30E+02  3.97E+02  8.19E+02 -1.91E+03 -2.56E+03
        -1.73E+03  8.75E+02  4.70E+02 -3.13E+03  9.61E+03
 
 OM44
+        6.18E+01 -1.63E+01 -1.82E+02  1.43E+02 -1.09E+02 -6.08E+01  1.17E+02 -1.08E+02  2.17E+02  5.16E+02 -3.18E+02 -3.19E+02
          2.13E+01  3.68E+02  9.43E+01 -2.69E+02  4.79E+02 -1.51E+02 -4.11E+02  2.09E+02  5.23E+02  9.10E+01 -3.50E+02  1.33E+02
         1.26E+02  3.00E+02  1.01E+01 -3.30E+02  3.51E+02  7.55E+02
 
 OM45
+       -1.08E+02 -8.93E+01  9.18E+01 -1.90E+02  1.27E+02  1.44E+02 -2.33E+02  2.75E+02 -1.61E+02 -9.23E+02  1.18E+02 -4.40E+02
         -6.32E+00 -8.99E+01 -1.91E+02  6.61E+02 -4.46E+02  7.89E+02 -2.96E+02 -6.07E+02 -9.97E+02 -2.80E+02  1.16E+03  2.98E+02
         5.05E+02  4.75E+01  6.44E+01  2.01E+02 -1.06E+03 -1.59E+02  1.59E+03
 
 OM46
+        6.98E+01  1.13E+02  1.35E+02 -1.16E+02  1.56E+02  3.76E+01 -1.42E+02  1.42E+01 -2.27E+02 -1.84E+02  1.26E+02  5.18E+02
         -1.87E+02 -5.94E+02  1.98E+02  5.04E+01 -3.68E+02  1.48E+02  6.42E+02 -1.13E+03 -1.25E+03  3.89E+02 -3.47E+02  3.75E+01
        -2.38E+02  3.06E+01  3.85E+02  1.94E+02 -4.50E+02 -5.88E+02  6.48E+02  1.73E+03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM47
+        1.45E+02  2.26E+02 -1.85E+02  2.41E+02 -2.51E+02 -1.62E+02  2.33E+02 -4.70E+02  1.87E+02  1.10E+03 -3.02E+02  5.99E+02
          7.11E+01  2.59E+02  6.79E+02 -7.07E+02  8.14E+02 -3.62E+02  1.15E+03  1.44E+02  4.64E+02  9.53E+02 -1.72E+03  2.40E+02
        -7.99E+02  2.00E+02  2.80E+02 -9.08E+02  1.34E+03  6.79E+02 -8.98E+02 -7.30E+02  2.64E+03
 
 OM48
+       -1.81E+02 -5.08E+02  1.90E+02 -2.05E+02  2.85E+02  1.14E+01 -4.42E+02  3.89E+02 -1.37E+02 -8.75E+02 -8.66E+01 -1.84E+03
          3.93E+02  1.82E+02 -1.13E+03  1.40E+03 -7.02E+02  7.47E+02 -2.34E+03  7.82E+02 -2.73E+02 -2.15E+03  1.70E+03 -6.32E+02
         2.20E+02 -1.09E+03 -3.92E+02  1.10E+03  1.12E+00 -4.55E+02  6.56E+02  3.77E+02 -1.82E+03  3.49E+03
 
 OM55
+       -2.01E+02 -1.23E+02  1.31E+02  6.15E+01  2.24E+02  1.39E+02 -1.65E+02  1.65E+02  1.46E+02  8.04E+02 -4.08E+02  9.98E+01
         -1.18E+03 -5.72E+02  7.19E+02 -1.01E+03  5.84E+02 -3.64E+02 -2.43E+02 -1.66E+03 -1.01E+03  5.24E+02 -1.51E+03 -2.12E+01
        -1.82E+02  1.57E+02 -3.24E+02 -4.61E+02  5.05E+02 -1.23E+02  2.01E+02  2.49E+02 -2.83E+02  2.67E+02  1.19E+03
 
 OM56
+       -3.21E+02 -4.08E+02  9.38E+01  1.53E+02  2.11E+02  1.64E+02 -7.35E+01  3.11E+02  1.00E+02  7.02E+02  5.63E+02 -1.49E+01
         -1.23E+03 -1.11E+03  3.64E+02 -5.07E+02  4.92E+02  6.77E+02 -9.77E+02 -2.44E+03 -1.14E+03  3.52E+02 -8.20E+02 -1.96E+02
        -1.68E+02 -1.26E+03 -1.43E+03  3.39E+02 -9.17E+02 -2.52E+01  9.17E+01  1.42E+02 -4.11E+02  5.88E+02  1.25E+03  4.15E+03
 
 OM57
+        1.05E+01  1.35E+02  1.02E+02 -2.30E+02 -3.51E+02 -1.10E+02 -1.29E+02 -2.12E+00 -6.09E+02 -1.81E+03  6.98E+02 -1.83E+02
          1.76E+03  3.77E+02 -1.27E+03  1.81E+03 -9.97E+02  1.30E+03 -2.21E+02  2.11E+03  4.68E+02 -2.08E+03  2.48E+03 -3.87E+02
         2.80E+02 -1.06E+03  1.24E+02  7.98E+02 -1.07E+03 -1.55E+02  2.98E+02 -2.41E+02  1.12E+01  4.38E+02 -7.83E+02 -1.01E+03
          2.79E+03
 
 OM58
+       -5.62E+01 -2.20E+02  3.40E+02  2.56E+02  2.74E+02  3.27E+02  8.91E+01 -4.59E+01  8.72E+02  2.62E+03 -1.02E+03  6.12E+02
         -3.62E+03 -9.43E+02  2.16E+03 -2.45E+03  1.38E+03 -1.84E+02  1.06E+03 -5.52E+03 -1.41E+03  2.92E+03 -3.71E+03 -3.98E+02
        -1.26E+03  2.62E+02 -8.70E+02 -8.52E+02  1.07E+03 -5.99E+01 -1.85E+02  4.81E+02  2.94E+02 -6.60E+02  1.66E+03  2.92E+03
         -2.61E+03  6.94E+03
 
 OM66
+       -1.36E+02 -2.04E+02  1.33E+01  3.10E+01  3.73E+01  7.88E+01 -9.17E+01  1.17E+02  7.92E+00  1.84E+02  1.67E+02 -1.85E+02
         -3.36E+02 -4.82E+02  1.27E+02 -1.18E+02  2.59E+02  8.74E+01 -5.62E+02 -2.87E+02 -5.78E+02 -2.03E+01 -3.34E+02 -3.37E+02
         9.61E+00 -8.17E+02 -6.63E+02  1.03E+02 -7.83E+01  4.92E+01 -5.46E+01 -3.12E+02 -1.47E+02  4.16E+02  5.75E+02  1.54E+03
         -1.23E+02  6.69E+02  1.10E+03
 
 OM67
+        1.49E+02  3.37E+02  1.64E+02 -1.54E+02 -7.37E+01 -2.03E+02  5.23E+01 -1.80E+02  2.48E+01 -4.85E+02 -5.15E+02  1.32E+02
          3.49E+02  5.29E+02 -9.36E+01  1.16E+02 -5.17E+02 -4.85E+02  2.37E+02  5.28E+02  7.84E+02 -5.10E+02  7.62E+02 -6.60E+01
         2.52E+02  3.30E+02  6.91E+01 -2.37E+02  2.80E+02 -1.76E+02 -1.03E+02  2.29E+02 -1.87E+02  7.24E+01 -3.28E+02 -7.02E+02
          5.20E+02 -6.99E+02 -5.14E+02  1.62E+03
 
 OM68
+       -1.93E+02 -7.30E+01 -1.03E+02  1.38E+01  2.39E+02  2.53E+02 -1.15E+02  1.70E+01 -8.67E+01  6.04E+02  5.96E+02  3.32E+02
         -7.33E+02 -1.34E+03  2.82E+02 -6.87E+02  7.56E+02  3.05E+02 -9.99E+01 -1.16E+03 -2.09E+03  7.50E+02 -1.76E+03 -5.44E+02
        -5.12E+02 -1.09E+03 -4.70E+02  3.66E+02 -1.74E+02 -1.77E+02  2.91E+02  4.25E+02 -1.13E+02  3.31E+02  9.59E+02  1.69E+03
         -4.62E+02  1.27E+03  9.54E+02 -9.45E+02  2.48E+03
 
 OM77
+        5.62E+01  8.00E+01 -1.27E+02  1.12E+02 -6.15E+01  3.14E+01  1.51E+02 -6.11E+01  1.82E+02  6.36E+02 -1.83E+02  5.14E+02
         -3.84E+02  1.50E+01  6.47E+02 -6.11E+02  4.54E+02 -4.29E+02  7.52E+02 -7.50E+02 -1.20E+02  1.19E+03 -1.05E+03  1.02E+02
        -6.45E+02  4.28E+02  1.56E+02 -8.54E+02  1.05E+03  1.81E+02 -2.04E+02 -6.20E+00  9.85E+02 -9.95E+02  5.70E+01 -7.48E+01
         -5.63E+02  7.34E+02 -9.43E+01 -3.45E+02  1.47E+02  9.43E+02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM78
+       -2.61E+02 -1.48E+02  4.74E+01 -4.31E+02  3.67E+01 -1.84E+02 -1.79E+02  3.25E+02 -7.20E+02 -2.44E+03  1.14E+03 -1.30E+03
          1.91E+03  3.04E+02 -2.50E+03  2.24E+03 -1.34E+03  8.05E+02 -2.54E+03  2.62E+03  8.53E+02 -3.44E+03  4.01E+03 -4.89E+02
         1.22E+03 -9.39E+02  7.76E+01  1.31E+03 -1.97E+03 -3.70E+02  7.38E+02  4.59E+01 -2.16E+03  2.42E+03 -5.38E+02 -3.52E+02
          1.88E+03 -3.37E+03  4.92E+01  1.03E+03 -6.13E+02 -1.69E+03  5.38E+03
 
 OM88
+        1.70E+02  4.22E+02 -4.08E+01  1.60E+02 -8.64E+01  1.31E+01  2.55E+02 -4.45E+02  4.10E+02  1.89E+03  6.10E+01  1.03E+03
         -9.71E+02 -4.86E+02  1.26E+03 -2.73E+03  1.36E+03 -1.02E+03  1.21E+03 -1.44E+03 -8.92E+02  1.78E+03 -3.74E+03  8.80E+02
        -1.93E+02  5.78E+02  3.09E+02 -7.88E+02 -5.35E+02  9.69E+01 -5.83E+02  3.15E+02  7.30E+02 -1.51E+03  6.10E+02  1.07E+02
         -1.31E+03  1.35E+03 -9.97E+01 -1.56E+02  8.27E+02  5.47E+02 -2.12E+03  3.07E+03
 
 SG11
+       -1.57E+03  6.76E+03 -1.75E+03 -1.21E+03 -2.96E+03  1.73E+03  9.45E+03 -2.88E+03  5.65E+03  1.78E+04  1.78E+04  6.97E+03
         -8.65E+03 -4.93E+03  1.05E+04 -1.92E+04  2.35E+04  1.26E+04 -1.93E+04 -1.55E+04  1.89E+04  7.03E+03 -9.50E+03 -9.87E+03
        -5.13E+02 -9.18E+03 -8.22E+03  4.44E+02 -1.50E+04  1.05E+04 -9.76E+03 -2.08E+04  6.29E+03 -1.34E+04 -6.85E+00  2.42E+04
         -1.64E+04  1.73E+04  1.10E+04 -2.73E+03  5.03E+03  6.65E+03 -5.51E+03  9.28E+03  2.42E+06
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG22
+       -2.72E+03 -3.30E+03  2.56E+03  4.37E+03  2.13E+03  8.52E+02 -1.53E+03  4.07E+03  6.35E+03  6.55E+03 -6.11E+03 -8.91E+03
         -6.17E+03  5.46E+03  7.34E+02 -2.39E+03  2.74E+02 -2.26E+03 -1.90E+04 -6.15E+03  9.78E+03 -5.44E+03  1.32E+04  2.39E+03
         2.14E+03  1.57E+04  4.22E+03 -7.92E+03 -5.91E+03  3.27E+03  6.12E+03 -3.11E+03 -6.85E+03  7.00E+03  3.75E+03  7.83E+02
         -9.79E+03  7.34E+03 -2.25E+03  3.07E+02 -1.73E+03 -9.30E+02 -1.41E+03  3.69E+01  1.84E+05  0.00E+00  8.40E+05
 
1
 
 
 #TBLN:      2
 #METH: MCMC Bayesian Analysis
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            3480
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      4
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     4
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example6hmto19_bayes.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 WISHART PRIOR DF INTERPRETATION (WISHTYPE):0
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 EM OR BAYESIAN METHOD USED:                MCMC BAYESIAN (BAYES)
 MU MODELING PATTERN (MUM):
 GRADIENT/GIBBS PATTERN (GRD):
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  0
 KEEP ITERATIONS (THIN):            1
 BURN-IN ITERATIONS (NBURN):                2000
 ITERATIONS (NITER):                        0
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       11456
 MC SAMPLES PER SUBJECT (ISAMPLE):          1
 RANDOM SAMPLING METHOD (RANMETHOD):        3U
 PROPOSAL DENSITY SCALING RANGE
              (ISCALE_MIN, ISCALE_MAX):     1.000000000000000E-06   ,1000000.00000000
 SAMPLE ACCEPTANCE RATE (IACCEPT):          0.400000000000000
 METROPOLIS HASTINGS SAMPLING FOR INDIVIDUAL ETAS:
 SAMPLES FOR GLOBAL SEARCH KERNEL (ISAMPLE_M1):          2
 SAMPLES FOR NEIGHBOR SEARCH KERNEL (ISAMPLE_M1A):       0
 SAMPLES FOR MASS/IMP/POST. MATRIX SEARCH (ISAMPLE_M1B): 2
 SAMPLES FOR LOCAL SEARCH KERNEL (ISAMPLE_M2):           2
 SAMPLES FOR LOCAL UNIVARIATE KERNEL (ISAMPLE_M3):       2
 PWR. WT. MASS/IMP/POST MATRIX ACCUM. FOR ETAS (IKAPPA): 1.00000000000000
 MASS/IMP./POST. MATRIX REFRESH SETTING (MASSREST):      1
 METROPOLIS HASTINGS POPULATION SAMPLING FOR NON-GIBBS
 SAMPLED THETAS AND SIGMAS:
 PROPOSAL DENSITY SCALING RANGE
              (PSCALE_MIN, PSCALE_MAX):   1.000000000000000E-02   ,1000.00000000000
 SAMPLE ACCEPTANCE RATE (PACCEPT):                       0.500000000000000
 SAMPLES FOR GLOBAL SEARCH KERNEL (PSAMPLE_M1):          1
 SAMPLES FOR LOCAL SEARCH KERNEL (PSAMPLE_M2):           -1
 SAMPLES FOR LOCAL UNIVARIATE KERNEL (PSAMPLE_M3):       1
 METROPOLIS HASTINGS POPULATION SAMPLING FOR NON-GIBBS
 SAMPLED OMEGAS:
 SAMPLE ACCEPTANCE RATE (OACCEPT):                       0.500000000000000
 SAMPLES FOR GLOBAL SEARCH KERNEL (OSAMPLE_M1):          -1
 SAMPLES FOR LOCAL SEARCH KERNEL (OSAMPLE_M2):           -1
 SAMPLES FOR LOCAL UNIVARIATE SEARCH KERNEL (OSAMPLE_M3):0
 DEG. FR. FOR T DIST.  PRIOR FOR THETAS (TTDF):        0.00000000000000
 DEG. FR. FOR LKJ CORRELATION PRIOR FOR OMEGAS (OLKJDF): 0.00000000000000
 WEIGHT FACTOR FOR STD PRIOR FOR OMEGAS (OVARF): 1.00000000000000
 DEG. FR. FOR LKJ CORRELATION PRIOR FOR SIGMAS (SLKJDF): 0.00000000000000
 WEIGHT FACTOR FOR STD PRIOR FOR SIGMAS (SVARF): 1.00000000000000

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8
 THETAS THAT ARE GIBBS SAMPLED:
   1   2   3   4   5   6   7   8
 THETAS THAT ARE METROPOLIS-HASTINGS SAMPLED:
 
 SIGMAS THAT ARE GIBBS SAMPLED:
   1   2
 SIGMAS THAT ARE METROPOLIS-HASTINGS SAMPLED:
 
 OMEGAS ARE GIBBS SAMPLED
 
 MONITORING OF SEARCH:

 Burn-in Mode
 iteration        -2000 MCMCOBJ=   -6804.69822710761     
 iteration        -1990 MCMCOBJ=   -6650.59904646967     
 iteration        -1980 MCMCOBJ=   -6669.35050846901     
 iteration        -1970 MCMCOBJ=   -6658.26908774069     
 iteration        -1960 MCMCOBJ=   -6632.97388907082     
 iteration        -1950 MCMCOBJ=   -6620.39942483948     
 iteration        -1940 MCMCOBJ=   -6608.10936939112     
 iteration        -1930 MCMCOBJ=   -6600.90692270886     
 iteration        -1920 MCMCOBJ=   -6648.64920730454     
 iteration        -1910 MCMCOBJ=   -6621.71926858169     
 iteration        -1900 MCMCOBJ=   -6553.58694466133     
 iteration        -1890 MCMCOBJ=   -6623.54764884824     
 iteration        -1880 MCMCOBJ=   -6580.22640140196     
 iteration        -1870 MCMCOBJ=   -6558.81360629037     
 iteration        -1860 MCMCOBJ=   -6620.33902605476     
 iteration        -1850 MCMCOBJ=   -6580.89258685854     
 iteration        -1840 MCMCOBJ=   -6605.75318602879     
 iteration        -1830 MCMCOBJ=   -6584.01639555086     
 iteration        -1820 MCMCOBJ=   -6531.63830436374     
 iteration        -1810 MCMCOBJ=   -6508.59509536754     
 iteration        -1800 MCMCOBJ=   -6562.40351998927     
 iteration        -1790 MCMCOBJ=   -6482.15238199588     
 iteration        -1780 MCMCOBJ=   -6555.84217502524     
 iteration        -1770 MCMCOBJ=   -6531.43449937059     
 iteration        -1760 MCMCOBJ=   -6555.22626010639     
 iteration        -1750 MCMCOBJ=   -6564.33498008596     
 iteration        -1740 MCMCOBJ=   -6564.40160619052     
 iteration        -1730 MCMCOBJ=   -6532.62470012466     
 iteration        -1720 MCMCOBJ=   -6544.74044738135     
 iteration        -1710 MCMCOBJ=   -6564.73488607305     
 iteration        -1700 MCMCOBJ=   -6602.35007987349     
 iteration        -1690 MCMCOBJ=   -6550.07817962979     
 iteration        -1680 MCMCOBJ=   -6536.31879940821     
 iteration        -1670 MCMCOBJ=   -6536.10171411276     
 iteration        -1660 MCMCOBJ=   -6549.67486039358     
 iteration        -1650 MCMCOBJ=   -6559.22599196363     
 iteration        -1640 MCMCOBJ=   -6530.27739518522     
 iteration        -1630 MCMCOBJ=   -6534.96305427894     
 iteration        -1620 MCMCOBJ=   -6563.84126048425     
 iteration        -1610 MCMCOBJ=   -6538.55584960950     
 iteration        -1600 MCMCOBJ=   -6548.09142936561     
 iteration        -1590 MCMCOBJ=   -6554.92590465467     
 iteration        -1580 MCMCOBJ=   -6537.53829707825     
 iteration        -1570 MCMCOBJ=   -6555.79905082050     
 iteration        -1560 MCMCOBJ=   -6491.72899695317     
 iteration        -1550 MCMCOBJ=   -6500.43441860961     
 iteration        -1540 MCMCOBJ=   -6596.75324176193     
 iteration        -1530 MCMCOBJ=   -6503.80577875503     
 iteration        -1520 MCMCOBJ=   -6558.48156477721     
 iteration        -1510 MCMCOBJ=   -6539.52861838156     
 iteration        -1500 MCMCOBJ=   -6533.35660974920     
 iteration        -1490 MCMCOBJ=   -6523.36004417712     
 iteration        -1480 MCMCOBJ=   -6544.29725644519     
 iteration        -1470 MCMCOBJ=   -6489.30939445032     
 iteration        -1460 MCMCOBJ=   -6532.84650489651     
 iteration        -1450 MCMCOBJ=   -6510.61501785137     
 iteration        -1440 MCMCOBJ=   -6545.06612349415     
 iteration        -1430 MCMCOBJ=   -6505.19240935793     
 iteration        -1420 MCMCOBJ=   -6546.21464418706     
 iteration        -1410 MCMCOBJ=   -6501.88124463855     
 iteration        -1400 MCMCOBJ=   -6557.94597158694     
 iteration        -1390 MCMCOBJ=   -6566.59028804776     
 iteration        -1380 MCMCOBJ=   -6505.93615105854     
 iteration        -1370 MCMCOBJ=   -6581.76368592821     
 iteration        -1360 MCMCOBJ=   -6490.66297299779     
 iteration        -1350 MCMCOBJ=   -6488.42592120995     
 iteration        -1340 MCMCOBJ=   -6495.06399192266     
 iteration        -1330 MCMCOBJ=   -6540.25565094557     
 iteration        -1320 MCMCOBJ=   -6489.02398721736     
 iteration        -1310 MCMCOBJ=   -6531.24267864620     
 iteration        -1300 MCMCOBJ=   -6543.05530338319     
 iteration        -1290 MCMCOBJ=   -6534.25661404501     
 iteration        -1280 MCMCOBJ=   -6442.32628679643     
 iteration        -1270 MCMCOBJ=   -6479.50035128108     
 iteration        -1260 MCMCOBJ=   -6484.46811026531     
 iteration        -1250 MCMCOBJ=   -6501.66577519103     
 iteration        -1240 MCMCOBJ=   -6493.26756077095     
 iteration        -1230 MCMCOBJ=   -6513.66179721493     
 iteration        -1220 MCMCOBJ=   -6481.87829978366     
 iteration        -1210 MCMCOBJ=   -6525.09639269064     
 iteration        -1200 MCMCOBJ=   -6454.71422199724     
 iteration        -1190 MCMCOBJ=   -6541.92504775091     
 iteration        -1180 MCMCOBJ=   -6518.59073757098     
 iteration        -1170 MCMCOBJ=   -6483.14694301589     
 iteration        -1160 MCMCOBJ=   -6533.47112999928     
 iteration        -1150 MCMCOBJ=   -6508.64220487172     
 iteration        -1140 MCMCOBJ=   -6530.53940664225     
 iteration        -1130 MCMCOBJ=   -6449.16575220687     
 iteration        -1120 MCMCOBJ=   -6557.93480148820     
 iteration        -1110 MCMCOBJ=   -6449.89793438356     
 iteration        -1100 MCMCOBJ=   -6468.64226771548     
 iteration        -1090 MCMCOBJ=   -6544.52328652071     
 iteration        -1080 MCMCOBJ=   -6483.69004843377     
 iteration        -1070 MCMCOBJ=   -6528.36714951391     
 iteration        -1060 MCMCOBJ=   -6491.84297851758     
 iteration        -1050 MCMCOBJ=   -6569.95768089925     
 iteration        -1040 MCMCOBJ=   -6545.81443168744     
 iteration        -1030 MCMCOBJ=   -6535.06279570608     
 iteration        -1020 MCMCOBJ=   -6480.09650057602     
 iteration        -1010 MCMCOBJ=   -6496.45183233563     
 iteration        -1000 MCMCOBJ=   -6478.90233910945     
 iteration         -990 MCMCOBJ=   -6500.41455036408     
 iteration         -980 MCMCOBJ=   -6475.38288688662     
 iteration         -970 MCMCOBJ=   -6522.61151666402     
 iteration         -960 MCMCOBJ=   -6497.70872092726     
 iteration         -950 MCMCOBJ=   -6499.03831975545     
 iteration         -940 MCMCOBJ=   -6501.19902753851     
 iteration         -930 MCMCOBJ=   -6528.30788174548     
 iteration         -920 MCMCOBJ=   -6457.33404032128     
 iteration         -910 MCMCOBJ=   -6526.63021107017     
 iteration         -900 MCMCOBJ=   -6444.83415095904     
 iteration         -890 MCMCOBJ=   -6481.76478451849     
 iteration         -880 MCMCOBJ=   -6539.32395179015     
 iteration         -870 MCMCOBJ=   -6502.43053677876     
 iteration         -860 MCMCOBJ=   -6480.77488127217     
 iteration         -850 MCMCOBJ=   -6470.66953846261     
 iteration         -840 MCMCOBJ=   -6465.15892232813     
 iteration         -830 MCMCOBJ=   -6438.04108723935     
 iteration         -820 MCMCOBJ=   -6499.15287050731     
 iteration         -810 MCMCOBJ=   -6527.02000006493     
 iteration         -800 MCMCOBJ=   -6453.63634508061     
 iteration         -790 MCMCOBJ=   -6516.34261722345     
 iteration         -780 MCMCOBJ=   -6522.50831839508     
 iteration         -770 MCMCOBJ=   -6497.21682795328     
 iteration         -760 MCMCOBJ=   -6536.42839259567     
 iteration         -750 MCMCOBJ=   -6552.93401993612     
 iteration         -740 MCMCOBJ=   -6478.91220360398     
 iteration         -730 MCMCOBJ=   -6439.79410379789     
 iteration         -720 MCMCOBJ=   -6534.62996640807     
 iteration         -710 MCMCOBJ=   -6495.96031280365     
 iteration         -700 MCMCOBJ=   -6492.35329421198     
 iteration         -690 MCMCOBJ=   -6499.59526368066     
 iteration         -680 MCMCOBJ=   -6497.32110952907     
 iteration         -670 MCMCOBJ=   -6496.26967057358     
 iteration         -660 MCMCOBJ=   -6471.29676265255     
 iteration         -650 MCMCOBJ=   -6457.51504820452     
 iteration         -640 MCMCOBJ=   -6434.72127559294     
 iteration         -630 MCMCOBJ=   -6535.40942123359     
 iteration         -620 MCMCOBJ=   -6508.63873951356     
 iteration         -610 MCMCOBJ=   -6450.91914288657     
 iteration         -600 MCMCOBJ=   -6512.09848212694     
 iteration         -590 MCMCOBJ=   -6462.48892847727     
 iteration         -580 MCMCOBJ=   -6504.02047141348     
 iteration         -570 MCMCOBJ=   -6513.43314385363     
 iteration         -560 MCMCOBJ=   -6473.67023575417     
 iteration         -550 MCMCOBJ=   -6517.90232992146     
 iteration         -540 MCMCOBJ=   -6447.56776409126     
 iteration         -530 MCMCOBJ=   -6412.16955091420     
 iteration         -520 MCMCOBJ=   -6511.49126714304     
 iteration         -510 MCMCOBJ=   -6503.10828237682     
 iteration         -500 MCMCOBJ=   -6550.59788166360     
 iteration         -490 MCMCOBJ=   -6491.85973766824     
 iteration         -480 MCMCOBJ=   -6480.01560457449     
 iteration         -470 MCMCOBJ=   -6466.53884082137     
 iteration         -460 MCMCOBJ=   -6427.34355882532     
 iteration         -450 MCMCOBJ=   -6494.28287294496     
 iteration         -440 MCMCOBJ=   -6501.47967241261     
 iteration         -430 MCMCOBJ=   -6495.08767641229     
 iteration         -420 MCMCOBJ=   -6432.18607161532     
 iteration         -410 MCMCOBJ=   -6490.63745177665     
 iteration         -400 MCMCOBJ=   -6500.01702985147     
 iteration         -390 MCMCOBJ=   -6511.83575633074     
 iteration         -380 MCMCOBJ=   -6508.49679886018     
 iteration         -370 MCMCOBJ=   -6468.30423604878     
 iteration         -360 MCMCOBJ=   -6451.11973388539     
 iteration         -350 MCMCOBJ=   -6458.78308681960     
 iteration         -340 MCMCOBJ=   -6463.71353268355     
 iteration         -330 MCMCOBJ=   -6478.06485285963     
 iteration         -320 MCMCOBJ=   -6476.51367964115     
 iteration         -310 MCMCOBJ=   -6528.23566830730     
 iteration         -300 MCMCOBJ=   -6489.70677907945     
 iteration         -290 MCMCOBJ=   -6490.54300369535     
 iteration         -280 MCMCOBJ=   -6404.93768867731     
 iteration         -270 MCMCOBJ=   -6443.21538743692     
 iteration         -260 MCMCOBJ=   -6482.76889207538     
 iteration         -250 MCMCOBJ=   -6509.05681731104     
 iteration         -240 MCMCOBJ=   -6516.23719510039     
 iteration         -230 MCMCOBJ=   -6468.83218589159     
 iteration         -220 MCMCOBJ=   -6463.77909721393     
 iteration         -210 MCMCOBJ=   -6508.19787156600     
 iteration         -200 MCMCOBJ=   -6516.35180624232     
 iteration         -190 MCMCOBJ=   -6446.10995037513     
 iteration         -180 MCMCOBJ=   -6467.42227382290     
 iteration         -170 MCMCOBJ=   -6500.41182146635     
 iteration         -160 MCMCOBJ=   -6501.92035962380     
 iteration         -150 MCMCOBJ=   -6473.08104205972     
 iteration         -140 MCMCOBJ=   -6459.45615199719     
 iteration         -130 MCMCOBJ=   -6455.98202335588     
 iteration         -120 MCMCOBJ=   -6478.53354224336     
 iteration         -110 MCMCOBJ=   -6491.69402596856     
 iteration         -100 MCMCOBJ=   -6494.11060143523     
 iteration          -90 MCMCOBJ=   -6464.67539158902     
 iteration          -80 MCMCOBJ=   -6447.37661646950     
 iteration          -70 MCMCOBJ=   -6511.31683153783     
 iteration          -60 MCMCOBJ=   -6528.46256476858     
 iteration          -50 MCMCOBJ=   -6438.22340232397     
 iteration          -40 MCMCOBJ=   -6524.26876053784     
 iteration          -30 MCMCOBJ=   -6524.97390296531     
 iteration          -20 MCMCOBJ=   -6515.68928045292     
 iteration          -10 MCMCOBJ=   -6465.16939006239     
 Sampling Mode
 iteration            0 MCMCOBJ=   -6479.79857755139     
 BURN-IN WAS NOT TESTED FOR CONVERGENCE
 
 #TERM:
 STATISTICAL PORTION WAS NOT PERFORMED
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         1568
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    2881.79124012985     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6479.79857755139     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -3598.00733742153     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           400
 NIND*NETA*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.150826563738     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6479.79857755139     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -5744.64775098765     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 PRIOR CONSTANT TO OBJECTIVE FUNCTION:    55.1779157436876     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6479.79857755139     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -6424.62066180770     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 #TERE:
 Elapsed estimation  time in seconds:   501.40
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 #OBJT:**************                       AVERAGE VALUE OF LIKELIHOOD FUNCTION                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -6479.799       **************************************************
 #OBJS:********************************************        0.000 (STD) **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8     
 
         4.05E+00 -2.27E+00  6.30E-01 -2.07E-01  2.35E+00  2.21E-01  3.73E+00 -6.40E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        3.57E-01
 
 ETA2
+       -3.69E-02  1.57E-01
 
 ETA3
+        5.57E-02 -2.66E-02  1.77E-01
 
 ETA4
+        6.40E-03  4.70E-02 -3.30E-02  2.43E-01
 
 ETA5
+        5.01E-02  1.95E-03  2.89E-02 -4.93E-02  1.88E-01
 
 ETA6
+       -7.20E-02  1.85E-03  1.32E-02  1.16E-02 -6.82E-02  2.25E-01
 
 ETA7
+        9.75E-02 -5.45E-02  8.39E-02 -6.80E-02  1.10E-01 -1.79E-02  3.38E-01
 
 ETA8
+        1.82E-01  3.72E-02  5.62E-02  1.24E-02  6.62E-02 -1.01E-01  1.62E-01  3.08E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1      EPS2     
 
 EPS1
+        9.69E-03
 
 EPS2
+        0.00E+00  2.43E-02
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        5.97E-01
 
 ETA2
+       -1.56E-01  3.96E-01
 
 ETA3
+        2.22E-01 -1.59E-01  4.21E-01
 
 ETA4
+        2.17E-02  2.41E-01 -1.59E-01  4.93E-01
 
 ETA5
+        1.93E-01  1.13E-02  1.58E-01 -2.31E-01  4.33E-01
 
 ETA6
+       -2.54E-01  9.86E-03  6.62E-02  4.96E-02 -3.32E-01  4.74E-01
 
 ETA7
+        2.81E-01 -2.37E-01  3.43E-01 -2.38E-01  4.35E-01 -6.51E-02  5.81E-01
 
 ETA8
+        5.50E-01  1.69E-01  2.40E-01  4.55E-02  2.75E-01 -3.85E-01  5.03E-01  5.55E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1      EPS2     
 
 EPS1
+        9.84E-02
 
 EPS2
+        0.00E+00  1.56E-01
 
1
 
 
 #TBLN:      3
 #METH: NUTS Bayesian Analysis
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            3480
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      4
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     4
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example6hmto19.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 WISHART PRIOR DF INTERPRETATION (WISHTYPE):0
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 EM OR BAYESIAN METHOD USED:                MCMC BAYESIAN (BAYES)
 MU MODELING PATTERN (MUM):
 GRADIENT/GIBBS PATTERN (GRD):
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  0
 KEEP ITERATIONS (THIN):            1
 BURN-IN ITERATIONS (NBURN):                250
 ITERATIONS (NITER):                        2000
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       11456
 MC SAMPLES PER SUBJECT (ISAMPLE):          1
 RANDOM SAMPLING METHOD (RANMETHOD):        3U
 PROPOSAL DENSITY SCALING RANGE
              (ISCALE_MIN, ISCALE_MAX):     1.000000000000000E-06   ,1000000.00000000
 SAMPLE ACCEPTANCE RATE (IACCEPT):          0.400000000000000
 METROPOLIS HASTINGS POPULATION SAMPLING FOR NON-GIBBS
 SAMPLED THETAS AND SIGMAS:
 PROPOSAL DENSITY SCALING RANGE
              (PSCALE_MIN, PSCALE_MAX):   1.000000000000000E-02   ,1000.00000000000
 SAMPLE ACCEPTANCE RATE (PACCEPT):                       0.500000000000000
 SAMPLES FOR GLOBAL SEARCH KERNEL (PSAMPLE_M1):          1
 SAMPLES FOR LOCAL SEARCH KERNEL (PSAMPLE_M2):           -1
 SAMPLES FOR LOCAL UNIVARIATE KERNEL (PSAMPLE_M3):       1
 METROPOLIS HASTINGS POPULATION SAMPLING FOR NON-GIBBS
 SAMPLED OMEGAS:
 SAMPLE ACCEPTANCE RATE (OACCEPT):                       0.500000000000000
 SAMPLES FOR GLOBAL SEARCH KERNEL (OSAMPLE_M1):          -1
 SAMPLES FOR LOCAL SEARCH KERNEL (OSAMPLE_M2):           -1
 SAMPLES FOR LOCAL UNIVARIATE SEARCH KERNEL (OSAMPLE_M3):0
 MASS/IMP./POST. MATRIX REFRESH SETTING (MASSREST):      0
 MASS MATRIX ACCUMULATION ITERATIONS (MADAPT):          200
 MASS MATRIX BLOCKING TYPE (NUTS_MASS):                 B
 MODEL PARAMETERS TRASNFORMED BY MASS MATRIX (NUTS_TRANSFORM=0)
 POWER TERM WEIGHTING FOR MASS MATRIX ACCUM. (KAPPA):   1.00000000000000
 NUTS SAMPLE ACCEPTANCE RATE (NUTS_DELTA):                   0.800000000000000
 NUTS GAMMA SETTING (NUTS_GAMMA):                            5.000000000000000E-02
 DEG. FR. FOR T DIST.  PRIOR FOR THETAS (TTDF):        0.00000000000000
 DEG. FR. FOR LKJ CORRELATION PRIOR FOR OMEGAS (OLKJDF): 8.00000000000000
 WEIGHT FACTOR FOR STD PRIOR FOR OMEGAS (OVARF): 1.00000000000000
 DEG. FR. FOR LKJ CORRELATION PRIOR FOR SIGMAS (SLKJDF): 0.00000000000000
 WEIGHT FACTOR FOR STD PRIOR FOR SIGMAS (SVARF): 1.00000000000000
 NUTS WARMUP METHOD (NUTS_TEST):       0
 NUTS MAXIMAL DEPTH SEARCH (NUTS_MAXDEPTH):       10
 NUTS STAGE I WARMUP ITERATIONS (NUTS_INIT):       7.500000000000000E-02
 NUTS STAGE II base WARMUP ITERATIONS (NUTS_BASE): 2.500000000000000E-02
 NUTS STAGE III FINAL ITERATIONS (NUTS_TERM): 5.000000000000000E-02
 INITIAL ITERATIONS FOR STEP NUTS SIZE ASSESSMENT (NUTS_STEPITER): 1
 INTERVAL ITERATIONS FOR STEP NUTS SIZE ASSESSMENT (NUTS_STEPINTER):0
 ETA PARAMETERIZATION (NUTS_EPARAM):0
 OMEGA PARAMETERIZATION (NUTS_OPARAM):1
 SIGMA PARAMETERIZATION (NUTS_SPARAM):1
 NUTS REGULARIZING METHOD (NUTS_REG): 0.00000000000000

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8
 THETAS THAT ARE GIBBS SAMPLED:
   1   2   3   4   5   6   7   8
 THETAS THAT ARE METROPOLIS-HASTINGS SAMPLED:
 
 SIGMAS THAT ARE GIBBS SAMPLED:
   1   2
 SIGMAS THAT ARE METROPOLIS-HASTINGS SAMPLED:
 
 OMEGAS ARE GIBBS SAMPLED
 
 MONITORING OF SEARCH:

 Burn-in Mode
 iteration         -250 MCMCOBJ=   -6485.77194769818     
 iteration         -245 MCMCOBJ=   -6602.54319673019     
 iteration         -240 MCMCOBJ=   -6603.75000689693     
 iteration         -235 MCMCOBJ=   -6594.28454485254     
 iteration         -230 MCMCOBJ=   -6584.66092781218     
 iteration         -225 MCMCOBJ=   -6516.22771080991     
 iteration         -220 MCMCOBJ=   -6581.90151074009     
 iteration         -215 MCMCOBJ=   -6558.98818572745     
 iteration         -210 MCMCOBJ=   -6647.00035661569     
 iteration         -205 MCMCOBJ=   -6563.94658281249     
 iteration         -200 MCMCOBJ=   -6560.74346092385     
 iteration         -195 MCMCOBJ=   -6623.80544710176     
 iteration         -190 MCMCOBJ=   -6626.48369500598     
 iteration         -185 MCMCOBJ=   -6623.11042108921     
 iteration         -180 MCMCOBJ=   -6565.60267962268     
 iteration         -175 MCMCOBJ=   -6575.58476079109     
 iteration         -170 MCMCOBJ=   -6560.97664776955     
 iteration         -165 MCMCOBJ=   -6593.88457281433     
 iteration         -160 MCMCOBJ=   -6622.89362644762     
 iteration         -155 MCMCOBJ=   -6604.60194947738     
 iteration         -150 MCMCOBJ=   -6608.94480442499     
 iteration         -145 MCMCOBJ=   -6607.36058937516     
 iteration         -140 MCMCOBJ=   -6605.86815653212     
 iteration         -135 MCMCOBJ=   -6551.13235279199     
 iteration         -130 MCMCOBJ=   -6560.92339817306     
 iteration         -125 MCMCOBJ=   -6590.63189534273     
 iteration         -120 MCMCOBJ=   -6636.61316095345     
 iteration         -115 MCMCOBJ=   -6604.58902326253     
 iteration         -110 MCMCOBJ=   -6622.94396917582     
 iteration         -105 MCMCOBJ=   -6608.11633705839     
 iteration         -100 MCMCOBJ=   -6535.33707014455     
 iteration          -95 MCMCOBJ=   -6637.60533517527     
 iteration          -90 MCMCOBJ=   -6602.02730840559     
 iteration          -85 MCMCOBJ=   -6603.20370034046     
 iteration          -80 MCMCOBJ=   -6548.09434306830     
 iteration          -75 MCMCOBJ=   -6598.23313065461     
 iteration          -70 MCMCOBJ=   -6635.13883375644     
 iteration          -65 MCMCOBJ=   -6575.85402064933     
 iteration          -60 MCMCOBJ=   -6619.91291480677     
 iteration          -55 MCMCOBJ=   -6576.58712339864     
 iteration          -50 MCMCOBJ=   -6582.22576100574     
 iteration          -45 MCMCOBJ=   -6556.15710515355     
 iteration          -40 MCMCOBJ=   -6606.06188760710     
 iteration          -35 MCMCOBJ=   -6626.25326480984     
 iteration          -30 MCMCOBJ=   -6579.77786288187     
 iteration          -25 MCMCOBJ=   -6541.43634294410     
 iteration          -20 MCMCOBJ=   -6562.51468357845     
 iteration          -15 MCMCOBJ=   -6571.72459256338     
 iteration          -10 MCMCOBJ=   -6582.38950855213     
 iteration           -5 MCMCOBJ=   -6637.49542145965     
 Sampling Mode
 iteration            0 MCMCOBJ=   -6605.42377830219     
 iteration            5 MCMCOBJ=   -6558.19371192720     
 iteration           10 MCMCOBJ=   -6642.66470386133     
 iteration           15 MCMCOBJ=   -6616.40091286170     
 iteration           20 MCMCOBJ=   -6587.47218960390     
 iteration           25 MCMCOBJ=   -6599.98838819967     
 iteration           30 MCMCOBJ=   -6655.48146399576     
 iteration           35 MCMCOBJ=   -6567.34608509322     
 iteration           40 MCMCOBJ=   -6609.23100695543     
 iteration           45 MCMCOBJ=   -6550.40223128854     
 iteration           50 MCMCOBJ=   -6657.98201216780     
 iteration           55 MCMCOBJ=   -6582.01495587338     
 iteration           60 MCMCOBJ=   -6554.70012487130     
 iteration           65 MCMCOBJ=   -6610.07307303271     
 iteration           70 MCMCOBJ=   -6600.37001658237     
 iteration           75 MCMCOBJ=   -6623.70527290689     
 iteration           80 MCMCOBJ=   -6578.21280843650     
 iteration           85 MCMCOBJ=   -6643.21012736548     
 iteration           90 MCMCOBJ=   -6619.64196745969     
 iteration           95 MCMCOBJ=   -6591.00302846871     
 iteration          100 MCMCOBJ=   -6521.42129514037     
 iteration          105 MCMCOBJ=   -6667.80622626404     
 iteration          110 MCMCOBJ=   -6625.70790134761     
 iteration          115 MCMCOBJ=   -6674.86158728055     
 iteration          120 MCMCOBJ=   -6587.83005168354     
 iteration          125 MCMCOBJ=   -6581.79807441908     
 iteration          130 MCMCOBJ=   -6658.39060389691     
 iteration          135 MCMCOBJ=   -6588.40375335902     
 iteration          140 MCMCOBJ=   -6608.98558341972     
 iteration          145 MCMCOBJ=   -6647.85805973316     
 iteration          150 MCMCOBJ=   -6634.07615488169     
 iteration          155 MCMCOBJ=   -6643.28903196509     
 iteration          160 MCMCOBJ=   -6569.56704748751     
 iteration          165 MCMCOBJ=   -6635.43301639685     
 iteration          170 MCMCOBJ=   -6633.63384914170     
 iteration          175 MCMCOBJ=   -6599.83895223263     
 iteration          180 MCMCOBJ=   -6588.30230288436     
 iteration          185 MCMCOBJ=   -6669.32224341400     
 iteration          190 MCMCOBJ=   -6587.37164207057     
 iteration          195 MCMCOBJ=   -6630.98610617270     
 iteration          200 MCMCOBJ=   -6613.25214718545     
 iteration          205 MCMCOBJ=   -6573.85793930073     
 iteration          210 MCMCOBJ=   -6558.20313132259     
 iteration          215 MCMCOBJ=   -6590.87744991424     
 iteration          220 MCMCOBJ=   -6652.30397299199     
 iteration          225 MCMCOBJ=   -6633.37714627749     
 iteration          230 MCMCOBJ=   -6597.84219675983     
 iteration          235 MCMCOBJ=   -6595.27732709214     
 iteration          240 MCMCOBJ=   -6637.71665890258     
 iteration          245 MCMCOBJ=   -6559.13526196059     
 iteration          250 MCMCOBJ=   -6565.14283196925     
 iteration          255 MCMCOBJ=   -6578.70778337095     
 iteration          260 MCMCOBJ=   -6565.67498475851     
 iteration          265 MCMCOBJ=   -6609.11060624427     
 iteration          270 MCMCOBJ=   -6492.27398144175     
 iteration          275 MCMCOBJ=   -6628.73834958619     
 iteration          280 MCMCOBJ=   -6613.55081155583     
 iteration          285 MCMCOBJ=   -6617.40385072083     
 iteration          290 MCMCOBJ=   -6553.49442724346     
 iteration          295 MCMCOBJ=   -6611.87490143744     
 iteration          300 MCMCOBJ=   -6651.90594468320     
 iteration          305 MCMCOBJ=   -6583.74252884745     
 iteration          310 MCMCOBJ=   -6643.18808796976     
 iteration          315 MCMCOBJ=   -6582.73487572303     
 iteration          320 MCMCOBJ=   -6617.11542645771     
 iteration          325 MCMCOBJ=   -6556.83746232411     
 iteration          330 MCMCOBJ=   -6581.90382054550     
 iteration          335 MCMCOBJ=   -6549.27296991843     
 iteration          340 MCMCOBJ=   -6600.67754757773     
 iteration          345 MCMCOBJ=   -6632.56793310986     
 iteration          350 MCMCOBJ=   -6614.25297750092     
 iteration          355 MCMCOBJ=   -6582.78863101701     
 iteration          360 MCMCOBJ=   -6572.28408591806     
 iteration          365 MCMCOBJ=   -6585.62427034562     
 iteration          370 MCMCOBJ=   -6574.66408145790     
 iteration          375 MCMCOBJ=   -6619.91975376898     
 iteration          380 MCMCOBJ=   -6652.78472093754     
 iteration          385 MCMCOBJ=   -6644.90967696470     
 iteration          390 MCMCOBJ=   -6634.00331467286     
 iteration          395 MCMCOBJ=   -6600.36848553454     
 iteration          400 MCMCOBJ=   -6589.04044208449     
 iteration          405 MCMCOBJ=   -6600.82916828781     
 iteration          410 MCMCOBJ=   -6624.06318659942     
 iteration          415 MCMCOBJ=   -6657.50927049130     
 iteration          420 MCMCOBJ=   -6624.16113088138     
 iteration          425 MCMCOBJ=   -6615.01607720281     
 iteration          430 MCMCOBJ=   -6667.02585719741     
 iteration          435 MCMCOBJ=   -6617.20614149527     
 iteration          440 MCMCOBJ=   -6650.58080548817     
 iteration          445 MCMCOBJ=   -6565.67410340741     
 iteration          450 MCMCOBJ=   -6666.20080658945     
 iteration          455 MCMCOBJ=   -6598.29750913880     
 iteration          460 MCMCOBJ=   -6587.26418793770     
 iteration          465 MCMCOBJ=   -6556.24373093262     
 iteration          470 MCMCOBJ=   -6625.11639042438     
 iteration          475 MCMCOBJ=   -6608.98250615103     
 iteration          480 MCMCOBJ=   -6599.59197715998     
 iteration          485 MCMCOBJ=   -6647.94547116062     
 iteration          490 MCMCOBJ=   -6663.39971005016     
 iteration          495 MCMCOBJ=   -6524.60945133437     
 iteration          500 MCMCOBJ=   -6619.26188187739     
 iteration          505 MCMCOBJ=   -6628.60412241015     
 iteration          510 MCMCOBJ=   -6623.80888655242     
 iteration          515 MCMCOBJ=   -6608.74425728514     
 iteration          520 MCMCOBJ=   -6643.45464443956     
 iteration          525 MCMCOBJ=   -6556.40478370022     
 iteration          530 MCMCOBJ=   -6561.31055137272     
 iteration          535 MCMCOBJ=   -6595.66444231912     
 iteration          540 MCMCOBJ=   -6629.80018318119     
 iteration          545 MCMCOBJ=   -6551.40559321548     
 iteration          550 MCMCOBJ=   -6626.24172382896     
 iteration          555 MCMCOBJ=   -6597.52130747443     
 iteration          560 MCMCOBJ=   -6614.18665607659     
 iteration          565 MCMCOBJ=   -6543.31328068385     
 iteration          570 MCMCOBJ=   -6578.16448152272     
 iteration          575 MCMCOBJ=   -6605.97585359608     
 iteration          580 MCMCOBJ=   -6557.49339557926     
 iteration          585 MCMCOBJ=   -6635.58862514757     
 iteration          590 MCMCOBJ=   -6601.88078525162     
 iteration          595 MCMCOBJ=   -6681.82178090930     
 iteration          600 MCMCOBJ=   -6580.68828896256     
 iteration          605 MCMCOBJ=   -6522.58889600686     
 iteration          610 MCMCOBJ=   -6578.04355579583     
 iteration          615 MCMCOBJ=   -6570.85561853335     
 iteration          620 MCMCOBJ=   -6569.34954971952     
 iteration          625 MCMCOBJ=   -6597.73848832546     
 iteration          630 MCMCOBJ=   -6555.99568405968     
 iteration          635 MCMCOBJ=   -6623.72379084425     
 iteration          640 MCMCOBJ=   -6637.16206174001     
 iteration          645 MCMCOBJ=   -6599.81670116449     
 iteration          650 MCMCOBJ=   -6643.78581480555     
 iteration          655 MCMCOBJ=   -6606.88687678624     
 iteration          660 MCMCOBJ=   -6584.12287959070     
 iteration          665 MCMCOBJ=   -6629.08089490500     
 iteration          670 MCMCOBJ=   -6586.89523583921     
 iteration          675 MCMCOBJ=   -6616.93447158692     
 iteration          680 MCMCOBJ=   -6571.20118631054     
 iteration          685 MCMCOBJ=   -6610.99218768154     
 iteration          690 MCMCOBJ=   -6587.11961260605     
 iteration          695 MCMCOBJ=   -6613.42468813636     
 iteration          700 MCMCOBJ=   -6580.12839521616     
 iteration          705 MCMCOBJ=   -6581.73430206883     
 iteration          710 MCMCOBJ=   -6652.98198536774     
 iteration          715 MCMCOBJ=   -6635.62282214522     
 iteration          720 MCMCOBJ=   -6678.74498702102     
 iteration          725 MCMCOBJ=   -6600.21013249745     
 iteration          730 MCMCOBJ=   -6627.97412585489     
 iteration          735 MCMCOBJ=   -6588.65546242090     
 iteration          740 MCMCOBJ=   -6604.25208505303     
 iteration          745 MCMCOBJ=   -6610.36905642966     
 iteration          750 MCMCOBJ=   -6543.05869537896     
 iteration          755 MCMCOBJ=   -6550.56382665999     
 iteration          760 MCMCOBJ=   -6598.10455267237     
 iteration          765 MCMCOBJ=   -6597.88122134838     
 iteration          770 MCMCOBJ=   -6585.68175568940     
 iteration          775 MCMCOBJ=   -6618.12262834407     
 iteration          780 MCMCOBJ=   -6615.47808407474     
 iteration          785 MCMCOBJ=   -6613.07904479968     
 iteration          790 MCMCOBJ=   -6638.91520063546     
 iteration          795 MCMCOBJ=   -6629.87833327326     
 iteration          800 MCMCOBJ=   -6638.40049722466     
 iteration          805 MCMCOBJ=   -6567.20079394396     
 iteration          810 MCMCOBJ=   -6598.06238151572     
 iteration          815 MCMCOBJ=   -6584.97090107210     
 iteration          820 MCMCOBJ=   -6567.05390563110     
 iteration          825 MCMCOBJ=   -6586.39822613779     
 iteration          830 MCMCOBJ=   -6594.40612203550     
 iteration          835 MCMCOBJ=   -6574.69401482675     
 iteration          840 MCMCOBJ=   -6646.84964505649     
 iteration          845 MCMCOBJ=   -6562.20656805634     
 iteration          850 MCMCOBJ=   -6649.68083674603     
 iteration          855 MCMCOBJ=   -6606.40131041774     
 iteration          860 MCMCOBJ=   -6686.78930634940     
 iteration          865 MCMCOBJ=   -6612.96433211318     
 iteration          870 MCMCOBJ=   -6579.17772675878     
 iteration          875 MCMCOBJ=   -6639.66753969329     
 iteration          880 MCMCOBJ=   -6594.57618634596     
 iteration          885 MCMCOBJ=   -6609.78640422631     
 iteration          890 MCMCOBJ=   -6642.21519423815     
 iteration          895 MCMCOBJ=   -6601.19030656071     
 iteration          900 MCMCOBJ=   -6584.07697480967     
 iteration          905 MCMCOBJ=   -6531.43511480334     
 iteration          910 MCMCOBJ=   -6596.53445303857     
 iteration          915 MCMCOBJ=   -6623.86545265321     
 iteration          920 MCMCOBJ=   -6588.91974175017     
 iteration          925 MCMCOBJ=   -6616.17866166304     
 iteration          930 MCMCOBJ=   -6566.23481077438     
 iteration          935 MCMCOBJ=   -6594.93318619951     
 iteration          940 MCMCOBJ=   -6571.18997256815     
 iteration          945 MCMCOBJ=   -6636.70566010126     
 iteration          950 MCMCOBJ=   -6577.40559969934     
 iteration          955 MCMCOBJ=   -6569.44375979824     
 iteration          960 MCMCOBJ=   -6615.07382892262     
 iteration          965 MCMCOBJ=   -6572.30624691006     
 iteration          970 MCMCOBJ=   -6566.54390738402     
 iteration          975 MCMCOBJ=   -6602.00838248484     
 iteration          980 MCMCOBJ=   -6638.53073397873     
 iteration          985 MCMCOBJ=   -6558.43722663064     
 iteration          990 MCMCOBJ=   -6632.06778447773     
 iteration          995 MCMCOBJ=   -6584.92010927623     
 Ending Mode
 iteration         1000 MCMCOBJ=   -6569.23447316088     
 BURN-IN WAS NOT TESTED FOR CONVERGENCE, AND WAS USER INTERRUPTED
 
 #TERM:
 STATISTICAL PORTION WAS NOT COMPLETED PRIOR TO USER INTERRUPT
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         1568
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    2881.79124012985     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6601.15400085214     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -3719.36276072229     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           400
 NIND*NETA*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.150826563738     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6601.15400085214     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -5866.00317428840     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 PRIOR CONSTANT TO OBJECTIVE FUNCTION:   -16.9020929079541     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6601.15400085214     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -6618.05609376009     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 #TERE:
 Elapsed estimation  time in seconds:  1478.68
 Elapsed covariance  time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              NUTS BAYESIAN ANALYSIS                            ********************
 #OBJT:**************                       AVERAGE VALUE OF LIKELIHOOD FUNCTION                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -6601.154       **************************************************
 #OBJS:********************************************       36.382 (STD) **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              NUTS BAYESIAN ANALYSIS                            ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8     
 
         3.91E+00 -2.21E+00  5.60E-01 -1.80E-01  2.27E+00  2.33E-01  3.71E+00 -7.04E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        2.67E-01
 
 ETA2
+       -3.01E-02  1.79E-01
 
 ETA3
+        2.79E-02 -1.17E-02  1.06E-01
 
 ETA4
+        2.33E-02  3.19E-02 -1.21E-02  2.51E-01
 
 ETA5
+        2.10E-02  1.92E-02 -1.38E-03 -2.29E-02  1.90E-01
 
 ETA6
+       -1.35E-02  1.15E-02  1.53E-02  1.15E-02 -5.60E-02  2.13E-01
 
 ETA7
+        1.04E-02 -3.61E-02  1.95E-02 -5.38E-02  1.77E-02  8.05E-03  2.25E-01
 
 ETA8
+        6.89E-02  6.40E-02  2.58E-02  3.42E-02 -6.48E-03 -4.08E-02  4.69E-02  1.92E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1      EPS2     
 
 EPS1
+        9.35E-03
 
 EPS2
+        0.00E+00  2.24E-02
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        5.14E-01
 
 ETA2
+       -1.37E-01  4.19E-01
 
 ETA3
+        1.70E-01 -8.56E-02  3.23E-01
 
 ETA4
+        9.09E-02  1.52E-01 -7.82E-02  4.98E-01
 
 ETA5
+        9.32E-02  1.07E-01 -9.39E-03 -1.06E-01  4.33E-01
 
 ETA6
+       -5.84E-02  6.14E-02  1.06E-01  4.86E-02 -2.84E-01  4.57E-01
 
 ETA7
+        4.58E-02 -1.78E-01  1.29E-01 -2.29E-01  8.60E-02  3.61E-02  4.72E-01
 
 ETA8
+        3.05E-01  3.48E-01  1.80E-01  1.55E-01 -3.34E-02 -2.04E-01  2.25E-01  4.36E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1      EPS2     
 
 EPS1
+        9.66E-02
 
 EPS2
+        0.00E+00  1.50E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              NUTS BAYESIAN ANALYSIS                            ********************
 ********************                STANDARD ERROR OF ESTIMATE (From Sample Variance)               ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8     
 
         7.21E-02  7.19E-02  5.14E-02  7.20E-02  6.61E-02  7.42E-02  6.29E-02  6.39E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        5.69E-02
 
 ETA2
+        3.06E-02  4.70E-02
 
 ETA3
+        2.16E-02  2.04E-02  2.94E-02
 
 ETA4
+        3.24E-02  2.90E-02  2.25E-02  5.36E-02
 
 ETA5
+        2.75E-02  2.36E-02  1.96E-02  2.82E-02  4.09E-02
 
 ETA6
+        2.99E-02  2.86E-02  2.02E-02  3.01E-02  2.73E-02  6.04E-02
 
 ETA7
+        3.01E-02  2.80E-02  2.01E-02  2.93E-02  2.49E-02  2.81E-02  4.83E-02
 
 ETA8
+        2.79E-02  2.51E-02  1.90E-02  2.78E-02  2.38E-02  2.72E-02  2.59E-02  3.90E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1      EPS2     
 
 EPS1
+        6.29E-04
 
 EPS2
+        0.00E+00  1.32E-03
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5      ETA6      ETA7      ETA8     
 
 ETA1
+        5.38E-02
 
 ETA2
+        1.33E-01  5.40E-02
 
 ETA3
+        1.25E-01  1.43E-01  4.37E-02
 
 ETA4
+        1.23E-01  1.31E-01  1.36E-01  5.25E-02
 
 ETA5
+        1.18E-01  1.27E-01  1.35E-01  1.24E-01  4.62E-02
 
 ETA6
+        1.20E-01  1.42E-01  1.33E-01  1.27E-01  1.27E-01  6.34E-02
 
 ETA7
+        1.21E-01  1.26E-01  1.28E-01  1.16E-01  1.17E-01  1.25E-01  4.94E-02
 
 ETA8
+        1.09E-01  1.13E-01  1.21E-01  1.19E-01  1.23E-01  1.25E-01  1.11E-01  4.34E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1      EPS2     
 
 EPS1
+        3.25E-03
 
 EPS2
+        0.00E+00  4.39E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              NUTS BAYESIAN ANALYSIS                            ********************
 ********************               COVARIANCE MATRIX OF ESTIMATE (From Sample Variance)             ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 TH 1
+        5.20E-03
 
 TH 2
+       -5.23E-04  5.17E-03
 
 TH 3
+        2.73E-04 -3.24E-04  2.64E-03
 
 TH 4
+        6.19E-04 -1.19E-04 -8.09E-05  5.19E-03
 
 TH 5
+       -3.54E-05  4.40E-04  2.20E-04 -5.67E-04  4.37E-03
 
 TH 6
+        6.76E-05 -1.06E-04 -1.81E-04  1.74E-04 -1.06E-03  5.51E-03
 
 TH 7
+        2.12E-04 -9.68E-04  4.99E-04 -1.02E-03  3.28E-04  5.26E-05  3.96E-03
 
 TH 8
+        1.36E-03  1.24E-03  5.35E-04  5.44E-04 -2.74E-04 -6.55E-04  1.05E-03  4.09E-03
 
 OM11
+       -1.45E-04  1.36E-05  5.89E-05 -1.63E-04 -9.81E-06  2.66E-04  3.55E-04 -1.35E-04  3.23E-03
 
 OM12
+       -2.83E-05  8.08E-06 -6.92E-05  4.95E-05  2.71E-06 -1.55E-04 -6.29E-05 -9.67E-05 -2.96E-04  9.36E-04
 
 OM13
+       -9.01E-05  3.97E-05 -5.05E-06  1.96E-05 -6.26E-05  5.82E-05 -1.13E-04 -7.14E-05  8.10E-05 -8.32E-05  4.66E-04
 
 OM14
+       -1.44E-04  7.54E-05  9.91E-05  6.58E-05 -3.87E-05  1.15E-04  2.93E-05  4.37E-05  5.46E-05  6.67E-05  3.08E-05  1.05E-03
 
 OM15
+       -9.39E-06  4.70E-05 -1.08E-04  1.01E-04  8.27E-05  8.36E-05  9.19E-05  8.64E-06  1.37E-04  1.97E-05 -2.72E-06 -2.46E-05
          7.56E-04
 
 OM16
+       -8.04E-05  7.46E-05  1.34E-05 -1.04E-04  1.04E-04  1.75E-05  9.90E-06  5.34E-05 -1.67E-04 -4.12E-05 -9.72E-06  8.19E-05
         -4.44E-05  8.95E-04
 
 OM17
+       -1.01E-04  1.01E-04 -1.18E-04 -4.11E-05  1.66E-04  9.67E-05 -1.97E-05  4.80E-05 -1.50E-04 -1.45E-04  3.36E-05 -1.61E-04
          7.61E-05  9.48E-05  9.04E-04
 
 OM18
+       -1.29E-04 -2.78E-05  4.52E-05 -4.99E-05  1.37E-05  6.62E-05  7.03E-05  4.73E-05  5.47E-04  1.00E-04  5.05E-05  1.44E-04
          4.93E-05 -1.34E-04  1.42E-04  7.81E-04
 
 OM22
+        3.40E-05 -5.75E-04  1.75E-04  9.51E-05  1.11E-04 -2.38E-05  1.02E-04  2.95E-05  1.15E-04 -3.01E-04  5.56E-06 -4.75E-05
          9.57E-05  6.03E-05  5.07E-05  3.72E-05  2.21E-03
 
 OM23
+       -1.11E-04  1.66E-04  2.70E-06 -1.40E-05  2.27E-05 -2.07E-05 -9.02E-06 -5.07E-05 -7.49E-05  7.26E-05 -4.09E-05  2.96E-05
          7.46E-06 -5.60E-05  8.39E-06  2.63E-06 -1.27E-04  4.18E-04
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM24
+        1.01E-04 -2.91E-04  1.04E-04 -4.71E-05  2.04E-05  6.02E-05  1.09E-04  3.05E-05 -1.13E-05 -4.83E-06  1.17E-05 -3.71E-05
         -1.77E-05 -2.83E-05 -5.28E-05 -8.51E-06  2.28E-04 -3.07E-05  8.42E-04
 
 OM25
+       -7.46E-05 -1.93E-05 -1.74E-05 -7.12E-05  4.39E-06 -8.65E-05  7.80E-05 -2.82E-06 -1.17E-04  6.40E-05 -5.78E-06  2.73E-05
         -4.50E-05 -1.06E-05  2.88E-06  2.99E-05  9.21E-05  4.95E-05 -4.51E-06  5.58E-04
 
 OM26
+       -1.21E-05  1.44E-04 -1.75E-05 -5.46E-05  4.73E-05 -1.20E-04  2.84E-05  5.78E-05  6.00E-05 -5.66E-05  2.20E-05  7.46E-05
         -2.09E-05 -6.03E-05 -5.44E-05 -2.48E-05  1.28E-05  5.79E-06  6.67E-05 -1.09E-04  8.20E-04
 
 OM27
+       -4.31E-05  4.49E-04 -4.49E-06 -6.62E-05 -8.40E-05  1.18E-04  4.71E-05  4.41E-05  4.14E-06  5.36E-05 -7.99E-06  1.92E-05
         -3.04E-05 -3.55E-05 -5.21E-05 -2.61E-05 -4.72E-04  1.26E-04 -1.77E-04  1.92E-05  4.18E-05  7.84E-04
 
 OM28
+       -3.38E-05 -6.32E-05  1.18E-04  5.85E-06  3.82E-05  1.55E-05  4.05E-05 -2.76E-05 -1.35E-05  1.56E-04 -7.14E-06 -8.16E-06
          3.93E-05 -2.69E-06 -3.88E-05  3.65E-06  4.56E-04  6.72E-05  1.04E-04 -7.78E-06 -8.97E-05  6.82E-05  6.28E-04
 
 OM33
+       -2.54E-05 -8.76E-07 -2.27E-04 -7.61E-05 -1.44E-04 -1.16E-04 -2.42E-05 -4.56E-05  2.09E-05 -1.85E-05  8.22E-05  4.05E-05
         -8.01E-06 -6.51E-05  5.03E-05  1.01E-04 -1.53E-05 -5.37E-05  3.54E-05  2.29E-05  6.72E-05 -1.32E-05 -3.32E-05  8.64E-04
 
 OM34
+       -1.92E-05  6.56E-05  1.02E-04  5.64E-05  9.07E-05  1.83E-05  1.17E-05 -2.17E-05 -5.44E-05  2.00E-05  6.69E-05  3.59E-05
          4.79E-06 -1.86E-05 -1.02E-05 -3.07E-05  2.22E-05  3.27E-05  4.95E-05  7.05E-06  1.40E-05  3.28E-05  3.98E-05  3.49E-06
         5.04E-04
 
 OM35
+        8.19E-06  4.19E-05 -1.94E-05 -4.18E-05  7.91E-05 -3.09E-05  8.50E-06  3.07E-05  4.14E-06  1.00E-05  9.23E-06 -5.33E-05
          2.36E-05 -7.18E-06  2.73E-05 -2.71E-06  4.40E-06  4.45E-05 -2.93E-05 -4.08E-06  3.54E-06  1.94E-05  6.81E-06 -4.75E-05
        -3.57E-05  3.85E-04
 
 OM36
+       -7.00E-05 -5.86E-05 -1.43E-05  2.69E-05 -6.75E-05  7.11E-05 -7.37E-05 -5.36E-05 -2.70E-05  9.68E-06 -4.19E-06  3.91E-05
          1.85E-05  1.82E-05 -2.08E-05  1.03E-05 -5.23E-05  1.93E-05 -1.58E-05  2.82E-05 -3.78E-05  1.07E-05 -1.16E-06  7.63E-06
         1.99E-05 -4.54E-05  4.10E-04
 
 OM37
+        5.39E-05 -1.06E-04 -1.12E-06  1.71E-05 -1.36E-04 -5.42E-05  3.74E-06 -2.57E-05  4.72E-06 -3.53E-05 -3.38E-05 -2.39E-06
         -2.60E-05 -1.26E-05  3.17E-05  4.08E-05 -5.43E-06 -8.85E-05 -1.55E-05 -1.76E-05  1.11E-05 -4.58E-05 -5.11E-05  7.64E-05
        -7.91E-05  4.13E-05  2.50E-05  4.05E-04
 
 OM38
+       -9.40E-05 -6.32E-05  7.04E-07 -8.54E-05 -2.82E-05 -9.62E-05  1.10E-05 -3.63E-05 -3.86E-05 -1.83E-05  1.25E-04  4.08E-05
         -1.06E-05 -7.17E-06  4.70E-05  6.40E-05  3.39E-05  7.92E-05  2.59E-05  2.68E-05  2.99E-05  1.05E-05  5.25E-05  1.55E-04
         8.53E-05 -2.90E-05 -4.73E-05  4.02E-05  3.59E-04
 
 OM44
+        9.05E-05 -1.03E-04  3.02E-04  2.60E-04  1.11E-04 -1.81E-04  1.97E-04  5.82E-05  8.62E-05  6.18E-05  8.71E-05  1.81E-04
          7.81E-05  7.73E-05 -7.17E-05 -3.12E-05  3.66E-05  1.08E-05  1.89E-04  5.82E-06  3.84E-05 -3.80E-05  4.86E-05 -1.50E-04
         8.70E-05  2.27E-05 -6.94E-05  2.27E-05  4.65E-05  2.87E-03
 
 OM45
+        4.19E-05 -7.66E-05 -5.75E-05  2.46E-05  1.48E-05  1.66E-05  1.27E-04 -1.11E-04 -1.42E-04  8.17E-05 -8.92E-06  2.35E-05
         -1.83E-06  5.08E-05 -3.15E-06 -3.17E-06  5.65E-05 -4.89E-06  3.91E-05  7.26E-05 -2.01E-05 -3.81E-05  7.52E-06  3.45E-05
        -1.35E-05 -9.45E-06 -3.72E-05 -1.16E-05  7.45E-07 -1.06E-04  7.93E-04
 
 OM46
+        4.99E-05  2.75E-05 -9.73E-05  5.49E-05  1.88E-06  2.62E-06  2.50E-05  1.54E-04  1.04E-04  3.44E-05  3.18E-05  3.93E-05
          2.65E-05  1.60E-04  1.58E-05  5.46E-06  2.99E-05 -5.80E-05  2.28E-05 -2.21E-05  7.23E-05 -7.93E-07  7.61E-06  3.74E-05
         3.42E-06 -1.24E-05  1.53E-05 -8.51E-06 -6.18E-06  1.86E-04 -4.39E-05  9.07E-04
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM47
+       -4.49E-05  9.74E-05 -6.91E-05 -6.47E-05 -4.29E-05  5.64E-05  1.03E-04  1.11E-05  2.96E-05  7.86E-05 -5.07E-05 -1.01E-05
          4.95E-05 -2.09E-05  9.96E-05  3.41E-05 -2.39E-05  4.76E-05 -1.48E-04  1.91E-06  1.14E-06  9.27E-05  2.06E-05 -8.97E-06
         4.52E-05  3.05E-05  1.30E-05 -4.54E-05 -3.07E-05 -2.77E-04  4.63E-05  4.66E-05  8.61E-04
 
 OM48
+        1.89E-06  1.16E-04  5.88E-05 -1.13E-04  8.09E-06  7.46E-05  1.22E-04  7.56E-05  1.23E-04  5.91E-05 -3.86E-06  2.14E-04
          5.25E-05 -2.04E-05 -5.85E-05  7.26E-05  1.81E-06  3.25E-05  2.36E-04 -3.17E-05  1.35E-05  3.46E-05  1.07E-04  2.62E-05
         8.81E-05 -4.42E-06  1.51E-06 -4.11E-05  2.00E-05  3.04E-04 -5.44E-05 -6.51E-05  1.26E-04  7.73E-04
 
 OM55
+       -7.73E-05  1.63E-04 -3.60E-05  1.21E-04  2.01E-05 -2.71E-04 -5.96E-05 -8.32E-05  1.98E-04 -6.19E-05  3.42E-05 -2.50E-06
          1.91E-04 -5.71E-05 -7.40E-05  2.21E-06  1.84E-05 -1.30E-05 -2.84E-05  2.81E-05  1.25E-05  4.22E-05 -4.02E-05  2.90E-05
         3.17E-05  1.12E-05  1.64E-05 -7.30E-05 -2.61E-05  7.09E-05 -8.24E-05 -1.53E-05  4.87E-06  9.48E-05  1.68E-03
 
 OM56
+       -7.49E-05 -3.42E-05 -1.67E-05 -5.91E-05  1.13E-05  8.43E-05 -2.23E-05 -3.07E-06 -2.21E-05 -1.95E-05 -7.09E-06 -8.71E-06
         -6.41E-06  8.44E-05  2.45E-05  2.11E-05  3.36E-05  1.01E-05  4.57E-05  2.33E-05  6.21E-05 -8.14E-06 -2.11E-05 -1.31E-05
         3.26E-05  9.40E-06 -1.06E-05  2.64E-06  2.03E-05  2.98E-05 -3.58E-07 -5.68E-05 -2.00E-05 -8.64E-06 -1.52E-04  7.46E-04
 
 OM57
+       -1.13E-05 -3.66E-05 -1.37E-05  7.61E-05 -2.51E-05  1.13E-05 -1.78E-05 -5.24E-05  7.92E-05 -2.21E-05  1.75E-05  2.83E-05
         -1.95E-06 -3.38E-05  1.06E-05  3.29E-05  1.59E-05  6.44E-06  7.50E-06 -8.34E-05  4.68E-05  4.49E-06  2.90E-05  9.04E-06
         1.27E-05 -8.59E-06  4.22E-06  4.95E-05  1.07E-06  6.68E-05 -1.11E-04 -1.40E-05 -1.04E-04 -2.05E-05  8.35E-05  4.54E-05
          6.18E-04
 
 OM58
+       -8.73E-05  7.24E-05 -5.10E-06  1.11E-04 -2.10E-05  4.88E-05  1.27E-05  3.44E-05 -3.49E-05  4.05E-05  1.44E-05 -2.19E-06
          1.10E-04 -3.33E-05  8.36E-05  5.44E-05  1.62E-05  3.48E-05 -2.72E-05  1.47E-04 -5.11E-05  1.94E-05  3.79E-05 -2.59E-05
         2.00E-06  6.43E-05  7.23E-06 -1.32E-05 -2.67E-06  6.48E-06  7.87E-05 -1.96E-06  2.49E-05 -7.23E-05 -1.19E-04 -9.21E-05
          8.27E-05  5.66E-04
 
 OM66
+        3.04E-05 -5.89E-05 -1.13E-04  2.27E-04 -3.33E-05  1.49E-04 -1.89E-04 -9.71E-06  2.20E-04  3.36E-05 -2.73E-05 -1.21E-05
         -4.74E-05  8.11E-05 -1.43E-04  1.52E-05  1.78E-04 -1.04E-04  6.87E-05 -3.86E-05  2.45E-05 -6.04E-05  2.19E-05  5.30E-05
         1.17E-05 -5.77E-05  9.88E-05  2.91E-05 -6.18E-06  3.24E-05  3.29E-05  2.23E-04 -3.16E-05  7.59E-06  2.30E-05 -4.49E-04
          9.32E-05 -2.19E-05  3.64E-03
 
 OM67
+       -3.73E-06  3.96E-05 -3.37E-05  9.94E-05  2.77E-06  1.94E-05 -4.54E-05  1.97E-06 -3.52E-05 -3.33E-05 -1.71E-05 -1.92E-05
         -6.05E-05 -3.04E-05 -1.27E-05 -1.50E-05 -8.02E-05 -1.09E-05 -4.84E-05 -2.01E-05 -1.70E-04 -4.84E-06 -2.39E-05 -1.66E-06
        -1.94E-05 -4.20E-05  2.49E-05  3.10E-05  5.62E-06 -4.07E-05 -1.18E-05 -9.55E-05  3.97E-05 -2.60E-05  2.55E-05  1.51E-05
         -6.62E-05 -9.23E-05  2.05E-04  7.87E-04
 
 OM68
+       -8.22E-05  7.24E-06  8.37E-06 -1.30E-04  1.83E-04  9.41E-06  5.73E-05 -3.77E-05 -3.78E-05 -8.71E-05 -8.77E-06  5.52E-05
         -4.15E-05  2.13E-04  2.21E-05 -5.29E-05  8.61E-06 -1.54E-05 -1.56E-05 -2.09E-05  1.79E-04 -1.59E-05 -1.21E-05  2.51E-05
        -3.15E-05  4.53E-06  5.96E-05  4.96E-06  7.65E-06 -2.29E-05  1.69E-05  1.10E-04 -1.35E-05 -4.22E-05 -3.63E-05  1.29E-06
         -5.65E-05 -1.08E-04 -3.58E-04  1.09E-04  7.40E-04
 
 OM77
+        1.23E-04 -1.92E-04  1.26E-04  7.08E-05  2.48E-05 -2.02E-04 -1.43E-04  1.16E-05  2.13E-05  1.29E-05 -2.73E-05 -1.25E-04
         -2.47E-05  4.56E-05 -1.34E-04 -2.19E-05  1.40E-04 -5.41E-05  6.25E-05  1.08E-05  8.14E-07 -3.36E-04  2.02E-05  1.82E-05
        -5.46E-05 -4.80E-05 -9.27E-06  7.05E-05 -8.10E-06 -1.17E-04  7.53E-05 -1.34E-05 -2.86E-04 -1.37E-04  1.49E-04 -5.33E-05
          1.35E-04  1.16E-05  2.03E-04  5.32E-05  4.17E-05  2.33E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM78
+        4.48E-05 -1.69E-04  5.74E-06 -8.90E-05  1.28E-06 -6.86E-06  3.24E-05  4.13E-05 -3.00E-05 -1.04E-05 -3.63E-06 -2.57E-05
          4.08E-05  6.65E-05  1.69E-04  6.18E-05 -1.13E-05 -1.98E-05  1.22E-06  1.05E-05  3.73E-05  9.15E-05 -3.41E-05  6.27E-05
        -6.08E-06 -1.81E-05  1.46E-05  5.22E-05  6.99E-05 -2.30E-05  3.03E-05  6.61E-06  6.89E-05 -5.83E-05 -4.29E-06  4.33E-05
         -8.42E-06  5.44E-05 -1.11E-04 -1.24E-04  1.86E-05  3.80E-04  6.71E-04
 
 OM88
+       -9.62E-05 -5.61E-05  7.70E-05 -2.08E-04 -4.22E-05 -1.48E-04  2.15E-04  2.14E-04  2.27E-04  1.03E-04  4.17E-05  5.12E-05
          5.02E-05 -5.64E-05  8.18E-05  3.45E-04  1.55E-04  5.61E-05  8.87E-05  4.07E-08 -1.29E-05  5.30E-05  3.83E-04  1.35E-04
         8.28E-05 -4.92E-05  4.82E-06 -1.10E-05  2.32E-04  1.78E-04 -3.09E-05 -3.97E-06  2.57E-05  2.63E-04 -1.57E-05  3.39E-05
          1.53E-05 -4.66E-06  5.68E-05 -7.37E-05 -1.78E-04  1.06E-04  3.23E-04  1.52E-03
 
 SG11
+        2.31E-06  1.72E-06  1.28E-06  3.24E-07  3.42E-06 -3.79E-06  9.17E-07  1.80E-06  9.17E-07  1.39E-06  1.19E-07  6.13E-07
          4.22E-07  2.24E-07 -1.02E-06 -5.63E-08 -2.38E-06 -2.30E-09 -8.51E-07 -2.39E-07 -9.27E-08  7.00E-07  1.34E-07  4.88E-08
        -9.14E-07  2.90E-07  4.70E-07 -6.79E-07  3.00E-07 -4.82E-07 -5.17E-07  2.32E-07  1.71E-06  3.23E-07  8.65E-07 -1.24E-06
          1.97E-08  2.03E-07  3.08E-06  1.16E-06 -2.76E-07 -5.26E-07 -5.06E-07  2.68E-07  3.96E-07
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG22
+        4.94E-06  5.47E-06  9.79E-07 -4.44E-06 -1.17E-07  4.38E-06  4.35E-06  2.54E-06  1.33E-06 -6.44E-07  4.73E-08  1.72E-06
         -6.27E-07 -5.50E-07  2.41E-06  2.73E-06 -4.02E-06 -6.87E-08 -2.92E-06  1.52E-06  2.87E-08  2.70E-06 -1.46E-06 -4.05E-09
        -6.55E-08 -1.40E-06  1.37E-08 -8.90E-07 -1.78E-06  2.46E-06 -2.83E-07  3.18E-06  3.33E-06 -1.90E-06 -2.09E-06 -7.55E-07
         -3.82E-07 -1.54E-06 -8.88E-06 -1.54E-06  3.44E-06 -8.05E-06 -2.32E-06 -3.82E-07  8.15E-09  0.00E+00  1.73E-06
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              NUTS BAYESIAN ANALYSIS                            ********************
 ********************              CORRELATION MATRIX OF ESTIMATE (From Sample Variance)             ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 TH 1
+        7.21E-02
 
 TH 2
+       -1.01E-01  7.19E-02
 
 TH 3
+        7.35E-02 -8.77E-02  5.14E-02
 
 TH 4
+        1.19E-01 -2.29E-02 -2.19E-02  7.20E-02
 
 TH 5
+       -7.42E-03  9.25E-02  6.48E-02 -1.19E-01  6.61E-02
 
 TH 6
+        1.26E-02 -1.99E-02 -4.75E-02  3.26E-02 -2.16E-01  7.42E-02
 
 TH 7
+        4.67E-02 -2.14E-01  1.54E-01 -2.25E-01  7.89E-02  1.13E-02  6.29E-02
 
 TH 8
+        2.95E-01  2.69E-01  1.63E-01  1.18E-01 -6.47E-02 -1.38E-01  2.60E-01  6.39E-02
 
 OM11
+       -3.53E-02  3.34E-03  2.01E-02 -3.99E-02 -2.61E-03  6.31E-02  9.93E-02 -3.72E-02  5.69E-02
 
 OM12
+       -1.28E-02  3.67E-03 -4.40E-02  2.25E-02  1.34E-03 -6.81E-02 -3.27E-02 -4.94E-02 -1.70E-01  3.06E-02
 
 OM13
+       -5.79E-02  2.56E-02 -4.55E-03  1.26E-02 -4.39E-02  3.63E-02 -8.32E-02 -5.17E-02  6.60E-02 -1.26E-01  2.16E-02
 
 OM14
+       -6.18E-02  3.24E-02  5.95E-02  2.82E-02 -1.81E-02  4.79E-02  1.44E-02  2.11E-02  2.97E-02  6.73E-02  4.40E-02  3.24E-02
 
 OM15
+       -4.73E-03  2.38E-02 -7.65E-02  5.09E-02  4.55E-02  4.09E-02  5.31E-02  4.91E-03  8.79E-02  2.35E-02 -4.59E-03 -2.76E-02
          2.75E-02
 
 OM16
+       -3.73E-02  3.47E-02  8.71E-03 -4.82E-02  5.26E-02  7.87E-03  5.26E-03  2.79E-02 -9.80E-02 -4.50E-02 -1.50E-02  8.45E-02
         -5.39E-02  2.99E-02
 
 OM17
+       -4.66E-02  4.67E-02 -7.66E-02 -1.90E-02  8.34E-02  4.33E-02 -1.04E-02  2.50E-02 -8.79E-02 -1.58E-01  5.18E-02 -1.65E-01
          9.20E-02  1.05E-01  3.01E-02
 
 OM18
+       -6.41E-02 -1.38E-02  3.15E-02 -2.48E-02  7.44E-03  3.19E-02  4.00E-02  2.65E-02  3.45E-01  1.17E-01  8.36E-02  1.59E-01
          6.42E-02 -1.60E-01  1.69E-01  2.79E-02
 
 OM22
+        1.00E-02 -1.70E-01  7.23E-02  2.81E-02  3.58E-02 -6.84E-03  3.47E-02  9.83E-03  4.31E-02 -2.09E-01  5.48E-03 -3.12E-02
          7.41E-02  4.29E-02  3.59E-02  2.84E-02  4.70E-02
 
 OM23
+       -7.54E-02  1.13E-01  2.57E-03 -9.54E-03  1.68E-02 -1.36E-02 -7.01E-03 -3.88E-02 -6.44E-02  1.16E-01 -9.26E-02  4.47E-02
          1.33E-02 -9.16E-02  1.36E-02  4.60E-03 -1.32E-01  2.04E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM24
+        4.81E-02 -1.40E-01  6.97E-02 -2.25E-02  1.06E-02  2.80E-02  5.95E-02  1.64E-02 -6.85E-03 -5.44E-03  1.88E-02 -3.95E-02
         -2.22E-02 -3.26E-02 -6.05E-02 -1.05E-02  1.68E-01 -5.17E-02  2.90E-02
 
 OM25
+       -4.37E-02 -1.13E-02 -1.43E-02 -4.19E-02  2.81E-03 -4.93E-02  5.25E-02 -1.87E-03 -8.69E-02  8.85E-02 -1.13E-02  3.56E-02
         -6.92E-02 -1.50E-02  4.06E-03  4.53E-02  8.30E-02  1.02E-01 -6.58E-03  2.36E-02
 
 OM26
+       -5.84E-03  7.01E-02 -1.19E-02 -2.65E-02  2.50E-02 -5.67E-02  1.58E-02  3.16E-02  3.69E-02 -6.46E-02  3.56E-02  8.04E-02
         -2.65E-02 -7.04E-02 -6.32E-02 -3.10E-02  9.54E-03  9.89E-03  8.03E-02 -1.61E-01  2.86E-02
 
 OM27
+       -2.13E-02  2.23E-01 -3.12E-03 -3.28E-02 -4.54E-02  5.68E-02  2.67E-02  2.46E-02  2.60E-03  6.26E-02 -1.32E-02  2.12E-02
         -3.95E-02 -4.24E-02 -6.19E-02 -3.34E-02 -3.59E-01  2.20E-01 -2.18E-01  2.91E-02  5.21E-02  2.80E-02
 
 OM28
+       -1.87E-02 -3.50E-02  9.19E-02  3.24E-03  2.31E-02  8.33E-03  2.57E-02 -1.72E-02 -9.49E-03  2.03E-01 -1.32E-02 -1.00E-02
          5.71E-02 -3.59E-03 -5.15E-02  5.22E-03  3.87E-01  1.31E-01  1.43E-01 -1.31E-02 -1.25E-01  9.72E-02  2.51E-02
 
 OM33
+       -1.20E-02 -4.14E-04 -1.50E-01 -3.60E-02 -7.40E-02 -5.33E-02 -1.31E-02 -2.42E-02  1.25E-02 -2.06E-02  1.30E-01  4.25E-02
         -9.91E-03 -7.40E-02  5.69E-02  1.23E-01 -1.11E-02 -8.93E-02  4.15E-02  3.30E-02  7.98E-02 -1.61E-02 -4.51E-02  2.94E-02
 
 OM34
+       -1.18E-02  4.06E-02  8.84E-02  3.49E-02  6.12E-02  1.10E-02  8.29E-03 -1.51E-02 -4.26E-02  2.91E-02  1.38E-01  4.94E-02
          7.76E-03 -2.77E-02 -1.51E-02 -4.90E-02  2.11E-02  7.12E-02  7.59E-02  1.33E-02  2.17E-02  5.22E-02  7.06E-02  5.29E-03
         2.25E-02
 
 OM35
+        5.79E-03  2.97E-02 -1.92E-02 -2.96E-02  6.10E-02 -2.12E-02  6.88E-03  2.45E-02  3.71E-03  1.67E-02  2.18E-02 -8.39E-02
          4.37E-02 -1.22E-02  4.62E-02 -4.95E-03  4.78E-03  1.11E-01 -5.14E-02 -8.81E-03  6.30E-03  3.54E-02  1.38E-02 -8.24E-02
        -8.09E-02  1.96E-02
 
 OM36
+       -4.79E-02 -4.02E-02 -1.37E-02  1.84E-02 -5.05E-02  4.73E-02 -5.79E-02 -4.14E-02 -2.34E-02  1.56E-02 -9.57E-03  5.96E-02
          3.33E-02  3.00E-02 -3.42E-02  1.82E-02 -5.50E-02  4.67E-02 -2.70E-02  5.89E-02 -6.52E-02  1.89E-02 -2.28E-03  1.28E-02
         4.39E-02 -1.14E-01  2.02E-02
 
 OM37
+        3.71E-02 -7.34E-02 -1.09E-03  1.18E-02 -1.02E-01 -3.63E-02  2.95E-03 -2.00E-02  4.12E-03 -5.73E-02 -7.77E-02 -3.67E-03
         -4.70E-02 -2.09E-02  5.24E-02  7.26E-02 -5.74E-03 -2.15E-01 -2.66E-02 -3.70E-02  1.92E-02 -8.13E-02 -1.01E-01  1.29E-01
        -1.75E-01  1.04E-01  6.13E-02  2.01E-02
 
 OM38
+       -6.88E-02 -4.63E-02  7.23E-04 -6.26E-02 -2.25E-02 -6.84E-02  9.20E-03 -2.99E-02 -3.58E-02 -3.16E-02  3.05E-01  6.65E-02
         -2.03E-02 -1.27E-02  8.25E-02  1.21E-01  3.81E-02  2.04E-01  4.71E-02  5.99E-02  5.51E-02  1.98E-02  1.11E-01  2.78E-01
         2.00E-01 -7.79E-02 -1.23E-01  1.05E-01  1.90E-02
 
 OM44
+        2.34E-02 -2.68E-02  1.10E-01  6.73E-02  3.14E-02 -4.56E-02  5.85E-02  1.70E-02  2.83E-02  3.77E-02  7.53E-02  1.04E-01
          5.30E-02  4.83E-02 -4.45E-02 -2.08E-02  1.45E-02  9.88E-03  1.22E-01  4.60E-03  2.50E-02 -2.54E-02  3.62E-02 -9.55E-02
         7.23E-02  2.16E-02 -6.40E-02  2.11E-02  4.59E-02  5.36E-02
 
 OM45
+        2.06E-02 -3.78E-02 -3.97E-02  1.22E-02  7.98E-03  7.92E-03  7.17E-02 -6.16E-02 -8.86E-02  9.48E-02 -1.47E-02  2.58E-02
         -2.36E-03  6.02E-02 -3.72E-03 -4.03E-03  4.27E-02 -8.50E-03  4.79E-02  1.09E-01 -2.49E-02 -4.83E-02  1.07E-02  4.17E-02
        -2.14E-02 -1.71E-02 -6.52E-02 -2.05E-02  1.40E-03 -7.01E-02  2.82E-02
 
 OM46
+        2.30E-02  1.27E-02 -6.28E-02  2.53E-02  9.46E-04  1.17E-03  1.32E-02  8.02E-02  6.07E-02  3.73E-02  4.89E-02  4.02E-02
          3.20E-02  1.77E-01  1.74E-02  6.49E-03  2.11E-02 -9.42E-02  2.61E-02 -3.10E-02  8.38E-02 -9.41E-04  1.01E-02  4.22E-02
         5.05E-03 -2.09E-02  2.50E-02 -1.40E-02 -1.08E-02  1.15E-01 -5.18E-02  3.01E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM47
+       -2.12E-02  4.61E-02 -4.58E-02 -3.06E-02 -2.21E-02  2.59E-02  5.61E-02  5.92E-03  1.78E-02  8.75E-02 -8.00E-02 -1.06E-02
          6.13E-02 -2.38E-02  1.13E-01  4.16E-02 -1.73E-02  7.93E-02 -1.74E-01  2.75E-03  1.35E-03  1.13E-01  2.81E-02 -1.04E-02
         6.86E-02  5.30E-02  2.18E-02 -7.69E-02 -5.53E-02 -1.76E-01  5.60E-02  5.27E-02  2.93E-02
 
 OM48
+        9.44E-04  5.82E-02  4.11E-02 -5.63E-02  4.40E-03  3.61E-02  6.95E-02  4.25E-02  7.77E-02  6.95E-02 -6.43E-03  2.38E-01
          6.86E-02 -2.45E-02 -7.00E-02  9.34E-02  1.39E-03  5.71E-02  2.93E-01 -4.82E-02  1.70E-02  4.45E-02  1.53E-01  3.20E-02
         1.41E-01 -8.11E-03  2.69E-03 -7.35E-02  3.80E-02  2.04E-01 -6.95E-02 -7.77E-02  1.54E-01  2.78E-02
 
 OM55
+       -2.62E-02  5.54E-02 -1.71E-02  4.12E-02  7.43E-03 -8.91E-02 -2.31E-02 -3.18E-02  8.51E-02 -4.94E-02  3.87E-02 -1.89E-03
          1.70E-01 -4.66E-02 -6.01E-02  1.93E-03  9.54E-03 -1.55E-02 -2.39E-02  2.91E-02  1.07E-02  3.69E-02 -3.92E-02  2.41E-02
         3.45E-02  1.40E-02  1.98E-02 -8.86E-02 -3.37E-02  3.23E-02 -7.15E-02 -1.24E-02  4.05E-03  8.33E-02  4.09E-02
 
 OM56
+       -3.80E-02 -1.74E-02 -1.19E-02 -3.00E-02  6.24E-03  4.16E-02 -1.30E-02 -1.76E-03 -1.42E-02 -2.34E-02 -1.20E-02 -9.84E-03
         -8.53E-03  1.03E-01  2.98E-02  2.76E-02  2.62E-02  1.80E-02  5.76E-02  3.62E-02  7.94E-02 -1.06E-02 -3.07E-02 -1.63E-02
         5.31E-02  1.75E-02 -1.91E-02  4.81E-03  3.92E-02  2.04E-02 -4.65E-04 -6.90E-02 -2.49E-02 -1.14E-02 -1.36E-01  2.73E-02
 
 OM57
+       -6.27E-03 -2.05E-02 -1.07E-02  4.25E-02 -1.53E-02  6.11E-03 -1.14E-02 -3.30E-02  5.60E-02 -2.91E-02  3.26E-02  3.51E-02
         -2.85E-03 -4.54E-02  1.42E-02  4.74E-02  1.36E-02  1.27E-02  1.04E-02 -1.42E-01  6.58E-02  6.45E-03  4.65E-02  1.24E-02
         2.28E-02 -1.76E-02  8.38E-03  9.90E-02  2.27E-03  5.02E-02 -1.58E-01 -1.86E-02 -1.42E-01 -2.96E-02  8.21E-02  6.69E-02
          2.49E-02
 
 OM58
+       -5.09E-02  4.23E-02 -4.17E-03  6.47E-02 -1.33E-02  2.76E-02  8.47E-03  2.26E-02 -2.58E-02  5.56E-02  2.81E-02 -2.84E-03
          1.68E-01 -4.68E-02  1.17E-01  8.18E-02  1.45E-02  7.15E-02 -3.95E-02  2.61E-01 -7.50E-02  2.92E-02  6.36E-02 -3.70E-02
         3.74E-03  1.38E-01  1.50E-02 -2.76E-02 -5.92E-03  5.08E-03  1.17E-01 -2.74E-03  3.56E-02 -1.09E-01 -1.22E-01 -1.42E-01
          1.40E-01  2.38E-02
 
 OM66
+        6.99E-03 -1.36E-02 -3.66E-02  5.22E-02 -8.34E-03  3.32E-02 -4.99E-02 -2.52E-03  6.40E-02  1.82E-02 -2.10E-02 -6.19E-03
         -2.85E-02  4.49E-02 -7.87E-02  9.00E-03  6.28E-02 -8.46E-02  3.92E-02 -2.71E-02  1.42E-02 -3.57E-02  1.45E-02  2.99E-02
         8.64E-03 -4.88E-02  8.09E-02  2.40E-02 -5.40E-03  1.00E-02  1.94E-02  1.23E-01 -1.78E-02  4.52E-03  9.32E-03 -2.72E-01
          6.21E-02 -1.52E-02  6.04E-02
 
 OM67
+       -1.84E-03  1.96E-02 -2.34E-02  4.92E-02  1.49E-03  9.34E-03 -2.57E-02  1.10E-03 -2.21E-02 -3.88E-02 -2.83E-02 -2.11E-02
         -7.84E-02 -3.62E-02 -1.51E-02 -1.91E-02 -6.09E-02 -1.90E-02 -5.94E-02 -3.04E-02 -2.12E-01 -6.16E-03 -3.40E-02 -2.02E-03
        -3.09E-02 -7.63E-02  4.39E-02  5.49E-02  1.06E-02 -2.71E-02 -1.49E-02 -1.13E-01  4.82E-02 -3.33E-02  2.22E-02  1.98E-02
         -9.49E-02 -1.38E-01  1.21E-01  2.81E-02
 
 OM68
+       -4.19E-02  3.70E-03  5.99E-03 -6.65E-02  1.02E-01  4.66E-03  3.35E-02 -2.17E-02 -2.44E-02 -1.05E-01 -1.49E-02  6.27E-02
         -5.55E-02  2.61E-01  2.70E-02 -6.97E-02  6.74E-03 -2.77E-02 -1.98E-02 -3.25E-02  2.30E-01 -2.08E-02 -1.78E-02  3.14E-02
        -5.16E-02  8.48E-03  1.08E-01  9.06E-03  1.48E-02 -1.57E-02  2.21E-02  1.34E-01 -1.69E-02 -5.58E-02 -3.26E-02  1.73E-03
         -8.35E-02 -1.67E-01 -2.18E-01  1.43E-01  2.72E-02
 
 OM77
+        3.52E-02 -5.53E-02  5.06E-02  2.04E-02  7.79E-03 -5.64E-02 -4.72E-02  3.76E-03  7.77E-03  8.74E-03 -2.62E-02 -8.02E-02
         -1.86E-02  3.16E-02 -9.22E-02 -1.63E-02  6.17E-02 -5.48E-02  4.46E-02  9.46E-03  5.89E-04 -2.49E-01  1.67E-02  1.28E-02
        -5.04E-02 -5.07E-02 -9.49E-03  7.25E-02 -8.85E-03 -4.54E-02  5.54E-02 -9.21E-03 -2.02E-01 -1.02E-01  7.56E-02 -4.05E-02
          1.13E-01  1.01E-02  6.98E-02  3.93E-02  3.18E-02  4.83E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM78
+        2.40E-02 -9.06E-02  4.31E-03 -4.77E-02  7.49E-04 -3.57E-03  1.99E-02  2.50E-02 -2.04E-02 -1.31E-02 -6.50E-03 -3.06E-02
          5.73E-02  8.59E-02  2.16E-01  8.53E-02 -9.29E-03 -3.74E-02  1.62E-03  1.72E-02  5.02E-02  1.26E-01 -5.25E-02  8.24E-02
        -1.05E-02 -3.55E-02  2.79E-02  1.00E-01  1.42E-01 -1.66E-02  4.16E-02  8.47E-03  9.07E-02 -8.10E-02 -4.05E-03  6.13E-02
         -1.31E-02  8.82E-02 -7.08E-02 -1.71E-01  2.64E-02  3.04E-01  2.59E-02
 
 OM88
+       -3.42E-02 -2.00E-02  3.84E-02 -7.40E-02 -1.64E-02 -5.10E-02  8.76E-02  8.58E-02  1.02E-01  8.63E-02  4.95E-02  4.05E-02
          4.68E-02 -4.83E-02  6.98E-02  3.17E-01  8.47E-02  7.04E-02  7.84E-02  4.42E-05 -1.16E-02  4.86E-02  3.92E-01  1.18E-01
         9.46E-02 -6.43E-02  6.11E-03 -1.40E-02  3.13E-01  8.50E-02 -2.81E-02 -3.38E-03  2.24E-02  2.42E-01 -9.82E-03  3.19E-02
          1.58E-02 -5.02E-03  2.41E-02 -6.74E-02 -1.68E-01  5.61E-02  3.20E-01  3.90E-02
 
 SG11
+        5.09E-02  3.81E-02  3.94E-02  7.16E-03  8.24E-02 -8.11E-02  2.32E-02  4.48E-02  2.56E-02  7.21E-02  8.78E-03  3.01E-02
          2.44E-02  1.19E-02 -5.38E-02 -3.20E-03 -8.04E-02 -1.79E-04 -4.66E-02 -1.61E-02 -5.15E-03  3.97E-02  8.50E-03  2.64E-03
        -6.47E-02  2.35E-02  3.69E-02 -5.36E-02  2.52E-02 -1.43E-02 -2.92E-02  1.22E-02  9.28E-02  1.84E-02  3.36E-02 -7.24E-02
          1.26E-03  1.36E-02  8.10E-02  6.59E-02 -1.61E-02 -1.73E-02 -3.10E-02  1.09E-02  6.29E-04
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG22
+        5.20E-02  5.77E-02  1.45E-02 -4.68E-02 -1.34E-03  4.48E-02  5.25E-02  3.02E-02  1.78E-02 -1.60E-02  1.66E-03  4.02E-02
         -1.73E-02 -1.40E-02  6.09E-02  7.41E-02 -6.50E-02 -2.55E-03 -7.64E-02  4.90E-02  7.62E-04  7.32E-02 -4.43E-02 -1.05E-04
        -2.22E-03 -5.43E-02  5.12E-04 -3.36E-02 -7.15E-02  3.49E-02 -7.62E-03  8.03E-02  8.61E-02 -5.17E-02 -3.87E-02 -2.10E-02
         -1.17E-02 -4.92E-02 -1.12E-01 -4.16E-02  9.59E-02 -1.27E-01 -6.81E-02 -7.44E-03  9.84E-03  0.00E+00  1.32E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              NUTS BAYESIAN ANALYSIS                            ********************
 ********************           INVERSE COVARIANCE MATRIX OF ESTIMATE (From Sample Variance)         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 TH 1
+        2.35E+02
 
 TH 2
+        5.63E+01  3.01E+02
 
 TH 3
+       -1.57E+00  3.95E+01  4.41E+02
 
 TH 4
+       -1.01E+01  3.76E+01  1.15E+01  2.31E+02
 
 TH 5
+       -2.03E+01 -5.31E+01 -2.29E+01  1.30E+01  2.77E+02
 
 TH 6
+       -2.09E+01 -2.36E+01  4.42E+00 -1.22E+01  6.24E+01  2.17E+02
 
 TH 7
+        3.18E+01  1.25E+02 -2.18E+01  8.26E+01 -4.58E+01 -3.37E+01  3.71E+02
 
 TH 8
+       -1.07E+02 -1.59E+02 -6.56E+01 -6.59E+01  6.80E+01  6.46E+01 -1.61E+02  4.17E+02
 
 OM11
+       -6.15E-01 -2.11E+01  2.88E+00 -7.39E-02  5.75E+00 -6.02E+00 -4.61E+01  4.26E+01  4.18E+02
 
 OM12
+        2.98E+00  1.69E-01  5.06E+01 -1.16E+01  7.72E+00  3.93E+01  6.47E+00  4.03E+01  2.00E+02  1.58E+03
 
 OM13
+        1.91E+01 -3.83E+01 -1.66E+01 -1.24E+01  4.44E+01 -2.55E+01  5.83E+01  4.29E+01 -3.64E+01  2.44E+02  2.75E+03
 
 OM14
+        3.67E+01 -3.52E-01 -3.41E+01 -2.29E+01 -6.25E+00 -3.38E+01  6.20E+00 -2.35E+01  3.70E+01  7.95E+00 -3.52E+01  1.21E+03
 
 OM15
+       -8.16E+00 -2.15E+01  7.11E+01 -2.81E+01 -3.00E+01 -2.61E+01 -5.79E+01  7.26E+00 -4.87E+01 -8.35E+01  3.06E+01  3.63E+01
          1.56E+03
 
 OM16
+        1.37E+01 -3.68E+01 -1.11E+01  2.16E+01  1.25E+00  6.96E+00 -1.21E+01 -5.70E+00  2.67E+01 -2.98E+01  8.50E+00 -1.74E+02
          1.05E+02  1.49E+03
 
 OM17
+        2.79E+01 -3.35E+01  4.69E+01 -6.76E+00 -5.97E+01 -3.31E+01  6.72E+00 -1.99E+01  1.51E+02  3.27E+02 -5.35E+01  3.03E+02
         -9.97E+01 -2.05E+02  1.54E+03
 
 OM18
+        2.76E+01  3.02E+01 -6.64E+01  3.33E+00 -3.36E+01 -2.71E+01  3.37E+01 -4.80E+01 -3.49E+02 -4.88E+02 -1.65E+02 -2.98E+02
          1.35E+01  3.06E+02 -4.48E+02  2.06E+03
 
 OM22
+        8.54E+00  3.36E+01 -1.01E+01 -9.73E+00 -4.70E+00  1.78E+00  6.42E+00 -2.22E+01  1.43E+01  4.14E+02  8.58E+01  2.30E+01
         -7.81E+01 -4.95E+01  6.46E+01 -1.62E+02  8.55E+02
 
 OM23
+        1.03E+01 -9.50E+01 -1.46E+01 -3.56E+01  2.63E+01  2.02E+01 -2.98E+01  8.50E+01  2.87E+01 -1.94E+01  6.65E+02 -1.01E+02
         -3.45E+01  1.61E+02 -1.35E+02 -2.81E+01  1.83E+02  3.38E+03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM24
+       -1.08E+00  5.52E+01 -3.06E+01  1.37E+01 -5.09E+00 -2.86E+01 -9.41E+00 -2.90E+01  2.21E+01  5.44E+01  5.36E+00  2.23E+02
          8.73E+01  6.67E+01  9.98E+01 -2.64E+01  2.04E+01 -7.80E+00  1.64E+03
 
 OM25
+        2.63E+00  1.82E+00  4.05E+01  3.20E+01 -2.57E+00  3.63E+01 -5.64E+01  6.66E+00  5.30E+01 -2.35E+02 -1.01E+01 -8.93E+01
          3.07E+02  1.23E+02 -1.75E+01  2.03E+01 -2.92E+02 -2.57E+02 -6.79E+01  2.39E+03
 
 OM26
+       -8.03E+00 -5.33E+01  1.27E+01  2.03E-01  3.19E+00  4.17E+01 -2.18E+01  6.11E+00 -3.54E+00 -6.94E+01 -8.68E+01 -1.48E+02
          7.41E+01  3.43E+02  2.86E+01  1.69E+02 -1.70E+02 -1.20E+02 -2.03E+02  4.31E+02  1.73E+03
 
 OM27
+       -8.06E+00 -1.41E+02 -3.49E+01 -7.46E+00  4.54E+01 -2.73E+01 -7.06E+01  3.29E+01  6.26E+00  2.60E+02  3.59E+01  1.33E+02
          9.40E+01  3.00E+00  3.09E+02 -8.40E+01  5.80E+02 -3.04E+02  4.52E+02 -2.82E+02 -2.75E+02  2.27E+03
 
 OM28
+       -1.31E+01  5.51E+00 -8.21E+01 -1.41E+01 -2.26E+01 -1.97E+01  2.02E+01  5.20E+01 -3.81E+01 -7.82E+02 -1.02E+02 -3.65E+01
         -2.66E+01  8.47E+01 -1.25E+02  5.24E+02 -9.33E+02 -2.96E+02 -3.65E+02  4.58E+02  5.58E+02 -1.00E+03  3.44E+03
 
 OM33
+       -1.58E+01 -2.47E+01  1.02E+02  4.91E+00  4.26E+01  3.76E+01 -1.60E+01  2.63E+01  1.33E+01  4.43E+00 -7.67E+01 -3.62E+01
          4.53E+00  1.20E+02 -5.13E+01 -8.61E+01  1.03E+01  2.60E+02 -5.28E+01 -6.46E+01 -5.51E+01 -9.84E+00  7.20E+01  1.41E+03
 
 OM34
+       -1.93E+01 -5.35E+01 -1.00E+02 -5.13E+01 -5.44E+01 -8.19E+00 -3.40E+01  5.20E+01  2.09E+01 -9.16E+01 -1.68E+02 -3.16E+01
          2.68E+01  8.06E+01 -1.86E+01  2.05E+02 -4.84E+01  1.20E+02 -8.47E+01  3.46E+01  1.40E+01 -9.28E+01  1.07E+02  4.62E+01
         2.37E+03
 
 OM35
+       -3.05E+00  3.98E+01  4.54E+01  4.51E+01 -6.33E+01 -3.85E+00  2.52E+01 -6.51E+01 -8.03E-01 -7.34E+01 -3.59E+02  1.35E+02
          1.93E+01  3.23E+01 -3.00E+01  1.33E+01 -6.37E+01 -5.91E+02  6.85E+01  1.75E+02  8.62E+01 -8.42E+01  3.20E+01  8.04E+01
         9.70E+01  3.06E+03
 
 OM36
+        3.55E+01  7.05E+01  7.07E+00  1.11E+01  2.56E+01 -1.81E+01  6.31E+01 -1.93E+01  5.05E+01  3.73E+01 -2.42E+02 -6.36E+01
         -1.11E+02  1.16E+01  1.04E+02 -4.06E+01  7.35E+01 -4.71E+02  1.54E+01 -1.32E+02  2.09E+02  4.87E+01 -1.75E+01 -6.97E+01
        -2.66E+02  4.65E+02  2.89E+03
 
 OM37
+       -3.94E+01 -5.31E+01 -5.47E+01 -3.54E+01  1.08E+02  5.29E+01 -5.29E+01  9.85E+01  1.71E+01  9.67E+01  5.79E+02 -2.46E+01
          2.93E+01  5.55E+01 -6.35E+01 -1.73E+02  6.71E+01  9.38E+02  7.72E+01 -5.04E+01 -1.06E+02  1.20E+02  3.40E+00 -1.25E+02
         4.62E+02 -6.75E+02 -4.86E+02  3.28E+03
 
 OM38
+        3.74E+01  9.62E+01  1.40E+01  6.30E+01  1.58E+01  4.15E+01  2.51E+01 -1.89E+01  9.62E+01  1.13E+02 -1.18E+03 -6.86E+01
         -2.38E+01 -5.80E+01 -2.14E+01 -7.85E+01 -2.85E+01 -1.33E+03 -5.66E+00 -1.31E+02 -1.68E+01  6.84E+01 -6.28E+01 -5.75E+02
        -6.85E+02  5.30E+02  8.68E+02 -9.35E+02  4.70E+03
 
 OM44
+       -6.09E+00 -2.18E+00 -3.97E+01 -2.94E+01 -3.45E-01  1.86E+01 -2.58E+01  2.41E+01 -1.15E+01 -4.88E+01 -8.38E+01 -2.51E+01
         -3.69E+01 -2.25E+01  6.07E+00  6.98E+01 -2.60E+00 -4.37E+01 -5.29E+00 -2.83E+01 -1.63E+01  4.55E+01  3.06E+01  8.58E+01
        -3.01E+01 -2.64E+01  8.19E+01 -6.72E+01  5.77E+00  4.28E+02
 
 OM45
+       -3.71E+01 -2.73E+01  2.99E+01 -3.00E+01  6.57E+00  9.93E+00 -9.07E+01  8.48E+01  4.95E+01 -1.19E+02 -5.23E+01 -8.51E+01
          1.37E+01 -7.22E+01  1.20E+01 -1.78E+01 -4.57E+01 -2.68E+00 -1.32E+02 -2.40E+01  4.84E+01 -1.12E+01  5.04E+01 -6.13E+01
         3.44E+01  8.19E+01  1.67E+02 -1.24E+00  5.05E+01  3.34E+01  1.45E+03
 
 OM46
+       -4.73E+00  1.53E+01  6.88E+01 -4.88E+00 -2.67E+00 -3.61E+00  9.83E+00 -6.33E+01 -5.14E+01 -1.09E+02 -8.20E+01 -6.51E+01
         -3.94E+01 -1.82E+02 -4.06E+01  1.30E+01 -1.69E+01  1.27E+02 -1.35E+02  5.74E+01 -6.84E+00 -1.01E+02  2.81E+01 -6.23E+01
        -1.87E+01  5.53E+01 -8.93E+00  5.63E+00  2.23E+01 -1.12E+02  1.23E+02  1.33E+03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM47
+        1.02E+01 -8.47E+00 -2.03E-01 -6.07E-01  3.96E+01  7.40E+00 -4.31E+01  2.06E+01 -2.54E+01 -1.04E+02  9.50E+01  1.08E+02
          6.53E-01  3.43E+01 -8.64E+01 -9.56E+00 -3.15E+01 -1.45E+02  3.85E+02  1.37E+01 -1.27E+02  1.35E+02 -1.38E+02  2.10E+01
        -1.78E+02 -9.72E+01  3.97E+00  7.29E+01  1.59E+02  1.86E+02 -1.06E+02 -1.92E+02  1.60E+03
 
 OM48
+       -2.05E+01 -4.21E+01  2.21E+01  3.69E+01 -1.52E+00 -2.58E+01 -7.01E+00 -1.08E+01 -2.91E+01 -2.91E+01  5.93E+01 -4.09E+02
         -1.20E+02 -3.67E+01 -3.22E+01 -1.16E+01  2.46E+01  2.15E+01 -6.56E+02  7.14E+01  9.24E+01 -2.04E+02  3.65E+01 -7.59E+01
        -1.54E+02 -4.38E+01 -2.92E+00  1.54E+01  1.26E+02 -1.98E+02  1.44E+02  2.75E+02 -5.00E+02  2.04E+03
 
 OM55
+        8.17E+00 -2.82E+01 -4.61E+00 -2.28E+01  1.49E+01  3.84E+01 -2.67E+00  3.33E+01 -2.08E+01  6.21E+01 -1.04E+01  1.90E+01
         -2.26E+02 -2.34E+01  4.35E+01 -3.14E+01 -1.41E+01  8.46E+01  1.42E+01 -1.90E+02 -4.91E+01 -5.86E+01  4.87E+01 -2.23E+01
        -3.07E+01 -1.15E+02 -5.55E+01  1.89E+02  1.46E+01 -1.14E+01  1.89E+01  4.03E-01 -1.23E+01 -6.63E+01  7.25E+02
 
 OM56
+        3.75E+01  1.55E+01  1.71E+01  5.64E+00 -1.40E+01 -4.09E+01  4.82E+01 -3.39E+01 -1.24E+01  4.43E+01  1.08E+02  6.93E+01
         -1.28E+02 -3.17E+02  6.52E+01 -1.43E+02 -3.39E+00  5.04E+01 -8.28E+01 -3.22E+02 -3.10E+02  3.89E+01 -9.59E+00  1.64E+01
        -1.31E+02 -1.79E+02 -9.24E+01  1.09E+02 -6.81E+01 -2.35E+01 -7.22E+01  4.44E+01  1.63E+01  7.23E+01  2.29E+02  1.78E+03
 
 OM57
+       -1.36E+01  2.15E+01  4.67E+01 -9.82E+00 -1.09E+00  6.80E+00 -3.66E+01  3.03E+01 -9.69E+00 -1.45E+01 -7.98E+01 -1.05E+02
          1.52E+02  6.91E+01 -1.11E+02 -1.73E+01 -3.77E+01 -1.55E+02 -1.30E+01  4.79E+02  4.40E+00 -1.34E+02 -5.66E+01 -2.91E+01
        -7.21E+01  2.05E+02  3.59E+01 -3.35E+02  1.29E+02 -2.33E+01  2.96E+02  7.64E+01  1.48E+02  7.30E+01 -1.85E+02 -3.34E+02
          2.02E+03
 
 OM58
+        5.01E+01 -4.59E+01 -2.39E+01 -4.45E+01  1.86E+01 -1.94E+01  3.37E+00 -2.18E+01  1.27E+01  1.08E+02  1.26E+01 -1.54E+01
         -4.24E+02 -6.04E+01 -9.34E+01 -2.06E+02  1.16E+02  8.44E+01 -9.23E+00 -8.30E+02 -1.22E+02  8.49E+01 -3.36E+02  5.04E+01
        -5.25E+01 -4.66E+02 -1.36E+02  2.30E+02 -2.69E+01 -3.62E+01 -2.48E+02 -2.62E+01 -6.77E+01  2.21E+02  3.18E+02  5.48E+02
         -5.66E+02  2.61E+03
 
 OM66
+        9.92E+00  8.65E+00  1.18E+01 -2.67E+00 -1.52E+01 -2.28E+01  2.30E+01 -1.03E+01 -2.29E+01  2.41E+00  7.50E+01  2.72E+01
         -3.05E-01 -1.35E+02  3.85E+01 -3.26E+01 -2.88E+01  8.46E+01 -2.04E+01 -6.03E+01 -1.46E+02 -1.56E+01  9.04E-01 -1.68E+01
        -2.14E+01 -3.54E+01 -1.42E+02  3.30E+01 -6.49E+01 -6.86E+00 -5.21E+01 -9.69E+01  5.42E+00  8.66E+00  4.61E+01  2.89E+02
         -9.99E+01  1.19E+02  3.85E+02
 
 OM67
+       -5.59E+00 -7.71E+00  3.40E+01 -2.63E+01 -5.86E+00  5.61E-01 -6.74E+00 -1.52E+01  2.44E+01  2.10E+01 -1.03E+01 -5.29E+01
          1.02E+02  2.24E+02 -8.68E+01  3.71E+01  1.08E+00  1.40E+01 -4.59E+01  1.90E+02  5.42E+02 -1.85E+02  1.92E+02  1.48E+01
         4.59E+01  2.35E+02  4.64E+01 -1.96E+02 -8.05E+01 -4.50E+01  6.89E+01  2.32E+02 -2.30E+02  1.55E+02 -7.17E+01 -2.65E+02
          2.28E+02  2.16E+01 -1.83E+02  1.73E+03
 
 OM68
+        3.75E+01  1.03E+01 -6.00E+00  1.44E+01 -6.65E+01 -2.82E+01 -1.15E+01  1.03E+00 -8.93E+00  2.06E+02  1.75E+02  5.35E+00
         -1.75E+01 -5.40E+02  4.49E+01 -1.53E+02  1.05E+02  1.28E+02  8.26E+01 -2.08E+02 -7.19E+02  1.84E+02 -4.99E+02 -6.76E+01
         6.42E+01 -2.40E+02 -4.60E+02  1.80E+02 -2.86E+02  8.61E+00 -1.01E+02 -2.38E+02  1.12E+02 -4.58E+01  1.03E+02  3.85E+02
         -2.89E+01  4.60E+02  3.19E+02 -5.82E+02  2.23E+03
 
 OM77
+       -2.22E+00 -2.39E+01 -3.32E+01 -2.44E+00  5.58E+00  5.77E+00  9.12E+00 -1.47E+00 -7.28E+00  2.21E+01  5.30E-01  1.06E+02
          2.12E+01 -2.07E+01  1.93E+02 -2.63E+01  7.78E+01 -9.33E+01  9.47E+01 -9.21E+01 -5.83E+01  4.38E+02 -2.20E+02 -1.20E+00
        -1.19E+01  2.77E+01  7.05E+01 -2.28E+01  1.15E+02  5.05E+01 -6.86E+01 -4.26E+01  2.34E+02 -4.95E+01 -4.69E+01  6.57E+01
         -1.57E+02  3.00E+01 -2.00E+01 -1.51E+02  6.35E+00  6.50E+02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      OM11      OM12      OM13      OM14  
             OM15      OM16      OM17      OM18      OM22      OM23      OM24      OM25      OM26      OM27      OM28      OM33  
            OM34      OM35      OM36      OM37      OM38      OM44      OM45      OM46      OM47      OM48      OM55      OM56  
             OM57      OM58      OM66      OM67      OM68      OM77      OM78      OM88      SG11      SG12      SG22  
 
 OM78
+       -2.98E+01  1.21E+02  9.84E+00  2.83E+01 -3.56E+01 -2.50E+01  5.95E+01 -3.48E+01  1.20E+00 -1.63E+02  5.13E+01 -1.52E+02
         -8.32E+01 -7.20E+01 -5.03E+02  1.99E+02 -2.09E+02  2.05E+02 -3.01E+02  2.00E+02  1.53E+02 -9.28E+02  9.06E+02 -2.77E+01
         9.09E+01  1.60E+02 -1.04E+02 -2.53E+02 -1.95E+02 -6.52E+01  3.19E+01  1.49E+02 -4.79E+02  4.49E+02 -3.85E+01 -1.60E+02
          2.18E+02 -2.77E+02  3.31E+01  5.14E+02 -3.10E+02 -5.91E+02  2.80E+03
 
 OM88
+        2.65E+01 -9.16E+00  2.25E+01  2.20E+01  1.96E+01  2.73E+01 -3.28E+01 -5.28E+01 -5.35E+00  1.42E+02  1.18E+02  1.01E+02
          1.39E+01 -4.24E+01  7.43E+01 -4.89E+02  1.74E+02  7.46E+01  1.55E+02 -7.54E+01 -1.74E+02  2.72E+02 -9.88E+02 -3.93E+01
        -5.72E+01  1.17E+00 -1.23E+02  1.67E+02 -5.29E+02 -3.28E+01 -1.13E+01 -4.62E+01  1.34E+02 -3.54E+02  2.97E+01  3.17E+01
         -1.32E+01  1.83E+02  1.07E+01 -1.13E+02  4.12E+02  7.49E+01 -7.62E+02  1.32E+03
 
 SG11
+       -1.40E+03 -6.86E+02 -1.72E+03 -5.29E+02 -2.03E+03  1.45E+03 -6.71E+02 -4.93E+01 -1.03E+03 -3.09E+03 -7.08E+02 -1.41E+03
         -1.58E+03 -1.83E+03  2.35E+03  1.82E+03  2.96E+03  4.92E+03  5.50E+02  4.05E+02 -1.00E+03  9.82E+01 -1.10E+03 -1.95E+02
         8.75E+03 -3.65E+03 -5.33E+03  7.33E+03 -9.34E+03  4.30E+02  2.00E+03  5.06E+02 -6.06E+03  3.05E+02 -4.58E+02  2.44E+03
         -1.75E+03 -1.20E+03 -1.63E+03 -3.88E+03  2.00E+03  8.95E+01  1.39E+03 -6.03E+01  2.73E+06
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG22
+       -8.56E+02 -6.22E+02 -2.71E+02  4.77E+02  1.86E+02 -4.18E+02 -3.95E+02  2.11E+02  2.81E+02  9.03E+02 -4.79E+02 -1.06E+03
         -1.61E+02  9.97E+02 -1.67E+03 -2.44E+03  5.90E+02  3.75E+02  2.18E+02 -2.84E+03  9.73E+02 -1.74E+03  1.43E+03 -3.54E+02
        -5.28E+02  2.89E+03  8.39E+02  4.49E+02  4.52E+03 -1.24E+03 -8.90E+01 -1.45E+03 -3.05E+03  3.84E+03  9.96E+02  1.75E+03
         -1.31E+03  3.36E+03  1.37E+03  2.27E+03 -3.20E+03  8.72E+02  3.99E+03 -1.82E+03 -2.09E+04  0.00E+00  6.50E+05
 
 Elapsed postprocess time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,     1987.328
Stop Time: 
Mon 10/31/2016 
12:52 AM
