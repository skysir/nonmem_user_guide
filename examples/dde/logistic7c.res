Wed 05/22/2019 
12:11 PM
;DDE
$PROBLEM LOGISTIC
; turn off second derivative assessments, sometimes even 1st derivatives if only simulating
$ABBR DERIV2=NO DERIV2=NOCOMMON 
$INPUT ID AMT TIME PRDV DV EVID MDV 
$DATA LOGISTIC6.csv IGNORE=C
$SUBROUTINES ADVAN16 TOL=6 ATOL=6 
$MODEL NCOMPARTMENTS=1

$PK
CALLFL=-2
MXSTEP=2000000000
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
KG=EXP(MU_1+ETA(1))
Y0=EXP(MU_2+ETA(2))
YSS=EXP(MU_3+ETA(3))
TAU1=EXP(MU_4+ETA(4))
; Initial conditions
A_0(1)=Y0



TSTOP=500.0
; INITIALIZING EQUATIONS FOR DDE COMPARTMENTS
TAU_1=1*TAU1

$DES
; AD_1_1 is the State value of A(1) delayed for time TAU1.
; AP_1_1 is the State value of A(1) in the past, for time delay TAU1.

; DELAY SETUP FOR EQUATION SET 1
 AP_1_1=Y0
; DELAY EQUATIONS FOR EQUATION SET 0 (BASE EQUATIONS)
;BASE EQUATIONS
 DADT(1)=KG*(1.0-AD_1_1/YSS)*A(1)

$ERROR
A1=A(1)


Y1=1.0
IPRED=A(1)
Y=IPRED*(1.0+EPS(1))

;$THETA
;-1.609     ; KG
;-0.0001     ; Y0
;2.3026    ; YSS
;1.609     ; TAU1

;$OMEGA (0.01)x4

;$SIGMA
;0.003


$THETA
-1.2
-2.34841E-02 
2.8
0.9

$OMEGA BLOCK(4)
0.1
6.83328E-04  0.1
1.05324E-03 -1.82494E-03  0.1
7.11753E-04  8.46056E-04  1.47772E-03  0.1

$SIGMA
0.01

$EST METHOD=ITS INTERACTION NOABORT SIGL=4 SIGLO=6 MCETA=10 NSIG=2 PRINT=1 NITER=100 CTYPE=3 FAST
$EST METHOD=IMP INTERACTION MAXEVAL=9999 NOABORT SIGL=6 NSIG=2 PRINT=1 NITER=100 CTYPE=3 MAPITER=0
$COV MATRIX=R UNCONDITIONAL
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  83) FUNCTIONS ARE USED IN ABBREVIATED CODE, BUT THE $SUBROUTINES
 RECORD DOES NOT INCLUDE THE "OTHER" OPTION.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       22 MAY 2019
Days until program expires :4025
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0 alpha version 7
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 LOGISTIC
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     2340
 NO. OF DATA ITEMS IN DATA SET:   7
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   5
 MDV DATA ITEM IS DATA ITEM NO.:  7
0INDICES PASSED TO SUBROUTINE PRED:
   6   3   2   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID AMT TIME PRDV DV EVID MDV
0FORMAT FOR DATA:
 (5E14.0/2E14.0)

 TOT. NO. OF OBS RECS:     2340
 TOT. NO. OF INDIVIDUALS:       30
0LENGTH OF THETA:   4
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS BLOCK FORM:
  1
  1  1
  1  1  1
  1  1  1  1
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
  -0.1200E+01 -0.2348E-01  0.2800E+01  0.9000E+00
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.1000E+00
                  0.6833E-03   0.1000E+00
                  0.1053E-02  -0.1825E-02   0.1000E+00
                  0.7118E-03   0.8461E-03   0.1478E-02   0.1000E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E-01
0COVARIANCE STEP OMITTED:        NO
 R MATRIX SUBSTITUTED:          YES
 S MATRIX SUBSTITUTED:           NO
 EIGENVLS. PRINTED:              NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:       FAST
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 TURN OFF Cholesky Transposition of R Matrix (CHOLROFF): NO
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):              -1
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0 alpha version 7

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF AND DELAY EQUATIONS (RADAR5, ADVAN16)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   4
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         COMP 1       ON         YES        YES        YES        YES
    2         OUTPUT       OFF        YES        NO         NO         NO
 INITIAL (BASE) TOLERANCE SETTINGS:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:   6
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          3
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   2

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE TIMES.
0PK SUBROUTINE INDICATES THAT COMPARTMENT AMOUNTS ARE INITIALIZED.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0DES SUBROUTINE USES FULL STORAGE MODE.
1
 
 
 #TBLN:      1
 #METH: Iterative Two Stage
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               FAST
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            528
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    10
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      6
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     4
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): logistic7c.ext
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
 CONVERGENCE TYPE (CTYPE):                  3
 CONVERGENCE INTERVAL (CINTERVAL):          1
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 ITERATIONS (NITER):                        100
 ANEAL SETTING (CONSTRAIN):                 1

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
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4
 THETAS THAT ARE SIGMA-LIKE:
 
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=   562.805606177590
 iteration            1 OBJ=   64.4389613526966
 iteration            2 OBJ=  -6.50600992166979
 iteration            3 OBJ=  -68.4469638621318
 iteration            4 OBJ=  -128.622588929444
 iteration            5 OBJ=  -188.233404439665
 iteration            6 OBJ=  -247.568135274584
 iteration            7 OBJ=  -306.699462303223
 iteration            8 OBJ=  -365.623645269133
 iteration            9 OBJ=  -424.302843936480
 iteration           10 OBJ=  -482.671632551876
 iteration           11 OBJ=  -540.631016037192
 iteration           12 OBJ=  -598.034594661951
 iteration           13 OBJ=  -654.652132624711
 iteration           14 OBJ=  -710.100646616067
 iteration           15 OBJ=  -763.685995491415
 iteration           16 OBJ=  -813.997495483805
 iteration           17 OBJ=  -857.712225869174
 iteration           18 OBJ=  -885.913013416366
 iteration           19 OBJ=  -888.637219382807
 iteration           20 OBJ=  -889.207310166459
 iteration           21 OBJ=  -889.294261328140
 iteration           22 OBJ=  -889.324649066779
 iteration           23 OBJ=  -889.333120641609
 iteration           24 OBJ=  -889.336447577829
 iteration           25 OBJ=  -889.337209980250
 iteration           26 OBJ=  -889.336775075725
 iteration           27 OBJ=  -889.336965723187
 Convergence achieved
 
 #TERM:
 OPTIMIZATION WAS COMPLETED


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -5.8543E-08  5.6712E-07  2.5035E-07 -6.5293E-07
 SE:             1.6571E-02  1.5077E-02  1.5074E-02  1.3000E-02
 N:                      30          30          30          30
 
 P VAL.:         1.0000E+00  9.9997E-01  9.9999E-01  9.9996E-01
 
 ETASHRINKSD(%)  1.3359E+00  4.1629E+00  3.7726E-01  9.8436E+00
 ETASHRINKVR(%)  2.6540E+00  8.1525E+00  7.5310E-01  1.8718E+01
 EBVSHRINKSD(%)  1.3351E+00  4.1624E+00  3.7643E-01  9.8376E+00
 EBVSHRINKVR(%)  2.6523E+00  8.1515E+00  7.5145E-01  1.8707E+01
 EPSSHRINKSD(%)  2.3561E+00
 EPSSHRINKVR(%)  4.6566E+00
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2340
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    4300.63233539787     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -889.336965723187     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       3411.29536967468     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           120
  
 #TERE:
 Elapsed estimation  time in seconds:    68.01
 Elapsed covariance  time in seconds:     0.08
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               ITERATIVE TWO STAGE                              ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -889.337       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               ITERATIVE TWO STAGE                              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4     
 
        -1.62E+00 -2.42E-02  2.29E+00  1.58E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        8.75E-03
 
 ETA2
+        6.62E-04  7.68E-03
 
 ETA3
+        1.05E-03 -1.80E-03  7.11E-03
 
 ETA4
+        5.56E-04  7.19E-04  1.49E-03  6.45E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        3.23E-03
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        9.36E-02
 
 ETA2
+        8.08E-02  8.76E-02
 
 ETA3
+        1.33E-01 -2.43E-01  8.43E-02
 
 ETA4
+        7.40E-02  1.02E-01  2.21E-01  8.03E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        5.68E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               ITERATIVE TWO STAGE                              ********************
 ********************                          STANDARD ERROR OF ESTIMATE (S)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4     
 
         2.89E-02  4.53E-02  2.67E-02  2.69E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        5.09E-03
 
 ETA2
+        2.50E-03  3.20E-03
 
 ETA3
+        3.59E-03  2.14E-03  3.27E-03
 
 ETA4
+        1.82E-03  4.52E-03  2.31E-03  2.86E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.34E-04
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        2.72E-02
 
 ETA2
+        3.10E-01  1.82E-02
 
 ETA3
+        4.58E-01  2.67E-01  1.94E-02
 
 ETA4
+        2.37E-01  6.38E-01  3.30E-01  1.78E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.18E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               ITERATIVE TWO STAGE                              ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (S)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      OM11      OM12      OM13      OM14      OM22      OM23      OM24      OM33  
             OM34      OM44      SG11  
 
 TH 1
+        8.37E-04
 
 TH 2
+       -5.09E-05  2.05E-03
 
 TH 3
+       -9.50E-05  1.83E-04  7.14E-04
 
 TH 4
+        1.30E-04  1.13E-04 -1.72E-04  7.25E-04
 
 OM11
+        6.50E-05  1.08E-04 -1.94E-05  3.02E-05  2.59E-05
 
 OM12
+       -4.15E-07 -5.45E-05 -1.61E-05  6.00E-07 -5.66E-06  6.24E-06
 
 OM13
+        2.30E-05 -8.20E-05 -4.92E-05  3.14E-05  3.96E-06  6.61E-07  1.29E-05
 
 OM14
+        6.29E-06  3.13E-06  8.60E-06  5.49E-06  1.69E-06 -7.41E-07  1.67E-06  3.31E-06
 
 OM22
+        3.05E-05 -5.84E-05 -3.65E-05  3.25E-05  5.71E-07  3.00E-06  6.40E-06  3.32E-07  1.02E-05
 
 OM23
+       -3.16E-05  1.59E-05  2.45E-05 -1.64E-05 -3.78E-06  1.68E-07 -2.92E-06 -2.39E-07 -3.87E-06  4.58E-06
 
 OM24
+        2.24E-05  1.56E-04  6.48E-06  2.88E-05  1.17E-05 -6.04E-06 -3.86E-06  1.75E-06 -1.92E-06  6.85E-07  2.04E-05
 
 OM33
+       -2.97E-05  6.93E-05  2.41E-06  1.76E-05  2.83E-06 -4.90E-07 -2.74E-06 -7.56E-07 -2.24E-06 -1.33E-07  1.23E-06  1.07E-05
 
 OM34
+        2.60E-05 -1.39E-05  2.11E-05 -2.10E-06  2.69E-06 -6.05E-07  6.64E-07  1.58E-06  5.42E-07 -1.60E-06 -2.03E-06 -9.65E-08
          5.32E-06
 
 OM44
+        1.51E-05  3.77E-05  2.31E-05 -3.23E-05  3.01E-06 -1.74E-06 -2.12E-06  8.34E-07 -8.03E-07  1.05E-07  4.00E-06 -7.92E-07
          2.41E-06  8.20E-06
 
 SG11
+        4.92E-07  1.16E-07 -3.21E-07 -4.84E-07  4.89E-08  2.47E-08 -4.07E-08 -1.42E-08 -2.37E-08 -1.98E-08 -2.38E-08 -7.94E-08
          2.49E-08  5.57E-08  1.81E-08
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               ITERATIVE TWO STAGE                              ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (S)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      OM11      OM12      OM13      OM14      OM22      OM23      OM24      OM33  
             OM34      OM44      SG11  
 
 TH 1
+        2.89E-02
 
 TH 2
+       -3.89E-02  4.53E-02
 
 TH 3
+       -1.23E-01  1.52E-01  2.67E-02
 
 TH 4
+        1.67E-01  9.29E-02 -2.39E-01  2.69E-02
 
 OM11
+        4.41E-01  4.71E-01 -1.43E-01  2.20E-01  5.09E-03
 
 OM12
+       -5.75E-03 -4.82E-01 -2.41E-01  8.92E-03 -4.45E-01  2.50E-03
 
 OM13
+        2.21E-01 -5.05E-01 -5.12E-01  3.25E-01  2.16E-01  7.37E-02  3.59E-03
 
 OM14
+        1.20E-01  3.80E-02  1.77E-01  1.12E-01  1.82E-01 -1.63E-01  2.56E-01  1.82E-03
 
 OM22
+        3.29E-01 -4.03E-01 -4.28E-01  3.77E-01  3.51E-02  3.76E-01  5.57E-01  5.71E-02  3.20E-03
 
 OM23
+       -5.10E-01  1.64E-01  4.28E-01 -2.85E-01 -3.47E-01  3.13E-02 -3.79E-01 -6.14E-02 -5.66E-01  2.14E-03
 
 OM24
+        1.72E-01  7.63E-01  5.37E-02  2.37E-01  5.10E-01 -5.35E-01 -2.38E-01  2.13E-01 -1.33E-01  7.09E-02  4.52E-03
 
 OM33
+       -3.15E-01  4.69E-01  2.76E-02  2.00E-01  1.70E-01 -6.00E-02 -2.33E-01 -1.27E-01 -2.14E-01 -1.91E-02  8.35E-02  3.27E-03
 
 OM34
+        3.89E-01 -1.33E-01  3.42E-01 -3.38E-02  2.29E-01 -1.05E-01  8.01E-02  3.77E-01  7.35E-02 -3.24E-01 -1.94E-01 -1.28E-02
          2.31E-03
 
 OM44
+        1.82E-01  2.91E-01  3.01E-01 -4.18E-01  2.06E-01 -2.44E-01 -2.06E-01  1.60E-01 -8.77E-02  1.72E-02  3.09E-01 -8.46E-02
          3.65E-01  2.86E-03
 
 SG11
+        1.26E-01  1.91E-02 -8.92E-02 -1.34E-01  7.15E-02  7.35E-02 -8.43E-02 -5.80E-02 -5.52E-02 -6.89E-02 -3.92E-02 -1.81E-01
          8.03E-02  1.45E-01  1.34E-04
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               ITERATIVE TWO STAGE                              ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (S)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      OM11      OM12      OM13      OM14      OM22      OM23      OM24      OM33  
             OM34      OM44      SG11  
 
 TH 1
+        3.29E+03
 
 TH 2
+       -2.63E+02  3.81E+03
 
 TH 3
+       -5.73E+02  1.15E+03  4.62E+03
 
 TH 4
+        7.17E+01 -7.20E+02 -9.59E+02  3.90E+03
 
 OM11
+       -4.02E+03 -9.03E+03 -3.06E+03  6.06E+03  1.37E+05
 
 OM12
+       -1.80E+04  1.20E+04  2.25E+04 -1.01E+04 -1.91E+04  5.71E+05
 
 OM13
+       -6.07E+03  2.25E+04  2.72E+04 -1.77E+04 -1.39E+05  2.73E+05  4.77E+05
 
 OM14
+        1.31E+04 -1.03E+04 -1.75E+04  1.11E+04  7.73E+04 -2.24E+05 -3.03E+05  6.48E+05
 
 OM22
+        6.02E+03  1.78E+02 -6.84E+03 -5.32E+03  3.20E+04 -2.44E+05 -1.58E+05  1.06E+05  3.49E+05
 
 OM23
+        2.23E+04 -1.36E+04 -3.28E+04  8.27E+03  7.17E+04 -3.83E+05 -2.65E+05  1.79E+05  3.05E+05  8.34E+05
 
 OM24
+       -1.00E+04 -1.03E+04  8.37E+03 -1.73E+04 -6.91E+04  2.15E+05  1.65E+05 -2.27E+05 -1.03E+05 -1.57E+05  3.62E+05
 
 OM33
+        1.34E+04 -1.64E+04 -3.88E+03 -4.83E+03 -1.35E+04 -1.22E+05 -4.89E+04  8.36E+04  6.37E+04  1.53E+05  5.03E+04  2.57E+05
 
 OM34
+       -1.65E+04  8.54E+03 -9.44E+03 -1.63E+04 -9.49E+04  1.63E+05  1.57E+05 -3.33E+05 -2.73E+04  2.57E+04  2.54E+05 -7.13E+04
          6.82E+05
 
 OM44
+        2.94E+03 -8.91E+03 -8.70E+03  2.35E+04  4.84E+04 -8.18E+04 -1.22E+05  1.20E+05 -3.81E+04  5.09E+04 -1.48E+05  1.81E+04
         -2.35E+05  3.38E+05
 
 SG11
+        2.80E+04 -5.04E+04  7.83E+04 -3.57E+04 -3.69E+05 -7.10E+05  5.35E+05  1.11E+05  6.35E+05  5.81E+05  6.06E+05  1.19E+06
         -2.84E+05 -5.75E+05  6.82E+07
 
1
 
 
 #TBLN:      2
 #METH: Importance Sampling
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               FAST
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    10
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      6
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     6
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): logistic7c.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 WISHART PRIOR DF INTERPRETATION (WISHTYPE):0
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 EM OR BAYESIAN METHOD USED:                IMPORTANCE SAMPLING (IMP)
 MU MODELING PATTERN (MUM):
 GRADIENT/GIBBS PATTERN (GRD):
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  3
 CONVERGENCE INTERVAL (CINTERVAL):          1
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 ITERATIONS (NITER):                        100
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       11456
 MC SAMPLES PER SUBJECT (ISAMPLE):          300
 RANDOM SAMPLING METHOD (RANMETHOD):        3U
 EXPECTATION ONLY (EONLY):                  0
 PROPOSAL DENSITY SCALING RANGE
              (ISCALE_MIN, ISCALE_MAX):     0.100000000000000       ,10.0000000000000
 SAMPLE ACCEPTANCE RATE (IACCEPT):          0.400000000000000
 LONG TAIL SAMPLE ACCEPT. RATE (IACCEPTL):   0.00000000000000
 T-DIST. PROPOSAL DENSITY (DF):             0
 NO. ITERATIONS FOR MAP (MAPITER):          0
 INTERVAL ITER. FOR MAP (MAPINTER):         0
 MAP COVARIANCE/MODE SETTING (MAPCOV):      1
 Gradient Quick Value (GRDQ):               0.00000000000000

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
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4
 THETAS THAT ARE SIGMA-LIKE:
 
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=  -889.347100101046 eff.=     300. Smpl.=     300. Fit.= 0.99093
 iteration            1 OBJ=  -889.537692450573 eff.=     124. Smpl.=     300. Fit.= 0.90143
 iteration            2 OBJ=  -888.835833968120 eff.=     119. Smpl.=     300. Fit.= 0.89973
 iteration            3 OBJ=  -889.841705252689 eff.=     121. Smpl.=     300. Fit.= 0.89894
 iteration            4 OBJ=  -889.344690763934 eff.=     122. Smpl.=     300. Fit.= 0.90101
 iteration            5 OBJ=  -889.089103311840 eff.=     120. Smpl.=     300. Fit.= 0.89961
 iteration            6 OBJ=  -889.038645183375 eff.=     117. Smpl.=     300. Fit.= 0.89850
 iteration            7 OBJ=  -889.778796865886 eff.=     122. Smpl.=     300. Fit.= 0.90001
 iteration            8 OBJ=  -889.427561013287 eff.=     122. Smpl.=     300. Fit.= 0.90171
 iteration            9 OBJ=  -889.548073742357 eff.=     120. Smpl.=     300. Fit.= 0.89830
 iteration           10 OBJ=  -889.651671983920 eff.=     121. Smpl.=     300. Fit.= 0.89968
 iteration           11 OBJ=  -889.368091809875 eff.=     121. Smpl.=     300. Fit.= 0.90085
 Convergence achieved
 iteration           11 OBJ=  -889.361865295049 eff.=     121. Smpl.=     300. Fit.= 0.90086
 
 #TERM:
 OPTIMIZATION WAS COMPLETED


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         6.6695E-05  1.8475E-04 -2.4982E-05 -6.7957E-04
 SE:             1.6453E-02  1.5132E-02  1.5064E-02  1.3250E-02
 N:                      30          30          30          30
 
 P VAL.:         9.9677E-01  9.9026E-01  9.9868E-01  9.5909E-01
 
 ETASHRINKSD(%)  1.5438E+00  3.8277E+00  4.7772E-01  9.5483E+00
 ETASHRINKVR(%)  3.0637E+00  7.5089E+00  9.5316E-01  1.8185E+01
 EBVSHRINKSD(%)  1.3643E+00  4.2092E+00  3.6953E-01  9.5577E+00
 EBVSHRINKVR(%)  2.7099E+00  8.2413E+00  7.3770E-01  1.8202E+01
 EPSSHRINKSD(%)  2.3522E+00
 EPSSHRINKVR(%)  4.6490E+00
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2340
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    4300.63233539787     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -889.361865295049     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       3411.27047010282     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           120
  
 #TERE:
 Elapsed estimation  time in seconds:    86.12
 Elapsed covariance  time in seconds:     7.33
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -889.362       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4     
 
        -1.62E+00 -2.37E-02  2.29E+00  1.57E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        8.67E-03
 
 ETA2
+        6.88E-04  7.68E-03
 
 ETA3
+        1.04E-03 -1.76E-03  7.11E-03
 
 ETA4
+        6.76E-04  7.48E-04  1.53E-03  6.66E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        3.23E-03
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        9.31E-02
 
 ETA2
+        8.43E-02  8.77E-02
 
 ETA3
+        1.32E-01 -2.38E-01  8.43E-02
 
 ETA4
+        8.89E-02  1.05E-01  2.23E-01  8.16E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        5.68E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                          STANDARD ERROR OF ESTIMATE (R)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4     
 
         1.73E-02  1.67E-02  1.55E-02  1.67E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        2.28E-03
 
 ETA2
+        1.58E-03  2.16E-03
 
 ETA3
+        1.46E-03  1.46E-03  1.84E-03
 
 ETA4
+        1.59E-03  1.53E-03  1.42E-03  2.06E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        9.71E-05
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.23E-02
 
 ETA2
+        1.92E-01  1.23E-02
 
 ETA3
+        1.81E-01  1.82E-01  1.09E-02
 
 ETA4
+        2.08E-01  2.08E-01  1.94E-01  1.26E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        8.54E-04
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (R)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      OM11      OM12      OM13      OM14      OM22      OM23      OM24      OM33  
             OM34      OM44      SG11  
 
 TH 1
+        2.98E-04
 
 TH 2
+        1.15E-05  2.80E-04
 
 TH 3
+        3.54E-05 -5.95E-05  2.39E-04
 
 TH 4
+        7.77E-06  3.84E-05  4.78E-05  2.79E-04
 
 OM11
+       -1.86E-07  5.57E-07 -1.15E-07  9.49E-07  5.21E-06
 
 OM12
+        3.99E-07 -8.76E-08 -1.30E-07 -1.76E-07  2.72E-07  2.48E-06
 
 OM13
+       -2.26E-07 -2.27E-08 -2.43E-07  5.77E-07  5.63E-07 -5.02E-07  2.13E-06
 
 OM14
+        4.36E-07 -5.58E-08  3.47E-07 -2.36E-06  9.56E-08  3.23E-07  3.95E-07  2.53E-06
 
 OM22
+       -3.06E-08 -4.65E-08  5.29E-09  1.46E-07 -2.19E-08  1.40E-07 -6.01E-08  7.23E-08  4.68E-06
 
 OM23
+       -1.35E-08 -9.90E-08  3.20E-07 -2.21E-07  3.11E-08  2.63E-07 -1.00E-08  5.99E-08 -9.52E-07  2.12E-06
 
 OM24
+        3.53E-07  2.86E-07 -1.31E-07 -8.76E-07 -3.23E-08  7.95E-08  6.29E-09  2.00E-07  6.63E-07  3.26E-07  2.33E-06
 
 OM33
+       -1.84E-07  3.36E-07 -1.33E-07  3.43E-07  1.08E-08 -1.58E-07  5.31E-07  1.10E-07  2.18E-07 -8.45E-07 -1.12E-07  3.38E-06
 
 OM34
+        4.04E-07 -2.76E-07  5.43E-07 -4.97E-07 -2.37E-08  4.30E-08  9.89E-08  2.76E-07 -1.41E-07  2.38E-07 -4.53E-07  6.77E-07
          2.03E-06
 
 OM44
+        1.04E-06 -5.61E-07  3.12E-07 -3.97E-06 -2.64E-08 -6.54E-09 -4.38E-08  2.82E-07  8.07E-08  1.38E-07  6.75E-07  1.16E-07
          8.76E-07  4.26E-06
 
 SG11
+       -1.59E-09 -5.45E-09 -8.81E-09  1.48E-09 -1.98E-10  1.67E-09 -5.77E-10  7.58E-10 -2.57E-09  1.57E-10 -1.33E-09 -1.18E-09
         -1.29E-10 -9.77E-10  9.42E-09
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (R)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      OM11      OM12      OM13      OM14      OM22      OM23      OM24      OM33  
             OM34      OM44      SG11  
 
 TH 1
+        1.73E-02
 
 TH 2
+        3.99E-02  1.67E-02
 
 TH 3
+        1.32E-01 -2.30E-01  1.55E-02
 
 TH 4
+        2.69E-02  1.37E-01  1.85E-01  1.67E-02
 
 OM11
+       -4.73E-03  1.46E-02 -3.26E-03  2.49E-02  2.28E-03
 
 OM12
+        1.47E-02 -3.32E-03 -5.32E-03 -6.66E-03  7.57E-02  1.58E-03
 
 OM13
+       -8.98E-03 -9.28E-04 -1.08E-02  2.36E-02  1.69E-01 -2.19E-01  1.46E-03
 
 OM14
+        1.59E-02 -2.10E-03  1.41E-02 -8.87E-02  2.64E-02  1.29E-01  1.70E-01  1.59E-03
 
 OM22
+       -8.19E-04 -1.28E-03  1.58E-04  4.04E-03 -4.44E-03  4.11E-02 -1.91E-02  2.10E-02  2.16E-03
 
 OM23
+       -5.37E-04 -4.06E-03  1.42E-02 -9.09E-03  9.36E-03  1.15E-01 -4.73E-03  2.59E-02 -3.02E-01  1.46E-03
 
 OM24
+        1.34E-02  1.12E-02 -5.53E-03 -3.43E-02 -9.27E-03  3.30E-02  2.82E-03  8.25E-02  2.00E-01  1.47E-01  1.53E-03
 
 OM33
+       -5.78E-03  1.09E-02 -4.69E-03  1.12E-02  2.57E-03 -5.45E-02  1.98E-01  3.75E-02  5.49E-02 -3.16E-01 -4.00E-02  1.84E-03
 
 OM34
+        1.64E-02 -1.16E-02  2.47E-02 -2.09E-02 -7.28E-03  1.91E-02  4.76E-02  1.22E-01 -4.56E-02  1.15E-01 -2.08E-01  2.58E-01
          1.42E-03
 
 OM44
+        2.91E-02 -1.62E-02  9.79E-03 -1.15E-01 -5.61E-03 -2.01E-03 -1.46E-02  8.62E-02  1.81E-02  4.59E-02  2.14E-01  3.05E-02
          2.98E-01  2.06E-03
 
 SG11
+       -9.46E-04 -3.35E-03 -5.87E-03  9.14E-04 -8.96E-04  1.09E-02 -4.07E-03  4.91E-03 -1.22E-02  1.11E-03 -8.96E-03 -6.64E-03
         -9.32E-04 -4.88E-03  9.71E-05
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (R)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      OM11      OM12      OM13      OM14      OM22      OM23      OM24      OM33  
             OM34      OM44      SG11  
 
 TH 1
+        3.43E+03
 
 TH 2
+       -2.69E+02  3.93E+03
 
 TH 3
+       -5.79E+02  1.17E+03  4.76E+03
 
 TH 4
+        2.52E+01 -7.42E+02 -9.85E+02  3.94E+03
 
 OM11
+        1.48E+02 -3.69E+02 -1.63E+01 -4.94E+02  2.01E+05
 
 OM12
+       -5.75E+02  5.00E+02  8.75E+02 -6.01E+02 -3.61E+04  4.52E+05
 
 OM13
+        9.27E+01  8.33E+02  1.28E+03 -1.77E+03 -6.48E+04  1.33E+05  5.61E+05
 
 OM14
+       -3.10E+02 -8.53E+02 -1.60E+03  3.86E+03  5.60E+03 -7.63E+04 -1.02E+05  4.36E+05
 
 OM22
+        2.18E+02  1.49E+02 -1.44E+02 -2.47E+02  9.81E+02 -2.55E+04 -1.96E+03  1.96E+03  2.56E+05
 
 OM23
+        6.22E+02 -1.57E+02 -7.93E+02  1.79E+02  3.56E+03 -7.17E+04 -5.21E+04  1.87E+04  1.46E+05  6.67E+05
 
 OM24
+       -5.16E+02 -5.78E+02  1.74E+02  5.73E+01  4.70E+03  4.73E+03  4.33E+03 -4.63E+04 -1.02E+05 -1.69E+05  5.59E+05
 
 OM33
+        4.11E+02 -4.76E+02  1.22E+01 -1.32E+02  7.94E+03 -1.41E+04 -9.44E+04  1.43E+04  2.21E+04  1.91E+05 -5.48E+04  3.87E+05
 
 OM34
+       -5.47E+02  2.44E+02 -6.81E+02 -8.11E+02  4.39E+03  2.48E+03  1.60E+04 -6.56E+04 -3.46E+04 -1.82E+05  2.14E+05 -1.70E+05
          6.94E+05
 
 OM44
+       -6.26E+02 -7.31E+01 -7.12E+02  3.52E+03 -2.39E+03  8.39E+03  1.10E+04 -6.55E+03  1.27E+04  3.30E+04 -1.21E+05  2.44E+04
         -1.62E+05  2.89E+05
 
 SG11
+       -3.09E+01  3.35E+03  5.15E+03 -2.01E+03  7.79E+03 -7.16E+04  9.93E+03 -3.63E+04  6.10E+04  3.72E+04  4.04E+04  3.96E+04
          5.84E+02  1.51E+04  1.06E+08
 
 Elapsed postprocess time in seconds:     0.00
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,      157.265
Stop Time: 
Wed 05/22/2019 
12:14 PM
