Fri 10/12/2018 
09:22 PM
$PROBLEM PDLIDR
$INPUT ID AMT TIME DV EVID  CMT
$DATA PDLIDR.csv IGNORE=C
$SUBROUTINES ADVAN16 TOL=4
$MODEL NCOMPARTMENTS=3 

$PK
MU_1=LOG(THETA(3))
MU_2=LOG(THETA(4))
MU_3=LOG(THETA(5))
MU_4=LOG(THETA(6))
MU_5=LOG(THETA(7))
MXSTEP=2000000000
KEL=THETA(1)
V=THETA(2)
K0=EXP(MU_1+ETA(1))
K1=EXP(MU_2+ETA(2))
SMAX=EXP(MU_3+ETA(3))
SC50=EXP(MU_4+ETA(4))
;  TAUy
TAU1=EXP(MU_5+ETA(5))
; Initial conditions
A_0(1)=0
A_0(2)=K0/K1
A_0(3)=K0*TAU1

$DES
; AD_x_y is the State value of A(x) delayed for time TAUy.  
; AP_x_y is the State value of A(x) in the past, for time delay TAUy.  
AP_2_1=K0/K1
;BASE EQUATIONS 
CC=A(1)/V
DADT(1)=-KEL*A(1)
DADT(2)=K0*(1+SMAX*CC/(SC50+CC))-K1*A(2)
DADT(3)=K1*A(2)-K1*AD_2_1

$ERROR

Y1=LOG(A(3))  
IF(CMT==3) IPRED=Y1
IF(CMT==3) Y=IPRED+EPS(1)

$THETA
0.25 FIX    ; 1: KEL
1    FIX    ; 2: V
(0,0.5,5)   ; 3: K0
(0,0.05,0.5); 4: K1
(0,50,500)  ; 5: SMAX
(0,1,10)    ; 6: SC50
(5,20,200)  ; 7: TR

$OMEGA BLOCK(5) VALUES(0.2,0.001)

$SIGMA
0.01
$EST METHOD=SAEM INTERACTION ISAMPLE=2 NBURN=200 NITER=200 PRINT=10 NOHABORT CTYPE=3 SIGL=3 RANMETHOD=3S2P
$EST METHOD=IMP INTERACTION MAPITER=0 NITER=200 ISAMPLE=300 CTYPE=3 PRINT=1
$COV UNCONDITIONAL MATRIX=R
$TABLE ID TIME IPRED EVID CMT NOPRINT ONEHEADER
FILE=PDLIDR.tab
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  3) THERE MAY BE AN ERROR IN THE ABBREVIATED CODE. THE FOLLOWING
 ONE OR MORE RANDOM VARIABLES ARE DEFINED WITH "IF" STATEMENTS THAT DO NOT
 PROVIDE DEFINITIONS FOR BOTH THE "THEN" AND "ELSE" CASES. IF ALL
 CONDITIONS FAIL, THE VALUES OF THESE VARIABLES WILL BE ZERO.
  
   IPRED Y

             
 (WARNING  83) FUNCTIONS ARE USED IN ABBREVIATED CODE, BUT THE $SUBROUTINES
 RECORD DOES NOT INCLUDE THE "OTHER" OPTION.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       12 OCT 2018
Days until program expires :4250
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0 alpha version 4
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 PDLIDR
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     3250
 NO. OF DATA ITEMS IN DATA SET:   7
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   4
 MDV DATA ITEM IS DATA ITEM NO.:  7
0INDICES PASSED TO SUBROUTINE PRED:
   5   3   2   0   0   0   6   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID AMT TIME DV EVID CMT MDV
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 IPRED
0FORMAT FOR DATA:
 (6E12.0,1F2.0)

 TOT. NO. OF OBS RECS:     3125
 TOT. NO. OF INDIVIDUALS:      125
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS BLOCK FORM:
  1
  1  1
  1  1  1
  1  1  1  1
  1  1  1  1  1
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.2500E+00     0.2500E+00     0.2500E+00
  0.1000E+01     0.1000E+01     0.1000E+01
  0.0000E+00     0.5000E+00     0.5000E+01
  0.0000E+00     0.5000E-01     0.5000E+00
  0.0000E+00     0.5000E+02     0.5000E+03
  0.0000E+00     0.1000E+01     0.1000E+02
  0.5000E+01     0.2000E+02     0.2000E+03
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.2000E+00
                  0.1000E-02   0.2000E+00
                  0.1000E-02   0.1000E-02   0.2000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.2000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.1000E-02   0.2000E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E-01
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
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING: (FPOSDEF):0
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    11456
 RANMETHOD:             3U
 MC SAMPLES (ESAMPLE):    300
 WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE
0-- TABLE   1 --
0RECORDS ONLY:    ALL
04 COLUMNS APPENDED:    YES
 PRINTED:                NO
 HEADERS:               ONE
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 IDFORMAT:
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID TIME IPRED EVID CMT
1DOUBLE PRECISION PREDPP VERSION 7.5.0 alpha version 4

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF AND DELAY EQUATIONS (RADAR5, ADVAN16)
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
    1            *           *           *           *           *
    2            *           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      5
   TIME DATA ITEM IS DATA ITEM NO.:          3
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   2
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:    6

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0PK SUBROUTINE INDICATES THAT COMPARTMENT AMOUNTS ARE INITIALIZED.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 
 
 #TBLN:      1
 #METH: Stochastic Approximation Expectation-Maximization
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            960
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
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      3
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     3
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): pdlidr.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 WISHART PRIOR DF INTERPRETATION (WISHTYPE):0
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 EM OR BAYESIAN METHOD USED:                STOCHASTIC APPROXIMATION EXPECTATION MAXIMIZATION (SAEM)
 MU MODELING PATTERN (MUM):
 GRADIENT/GIBBS PATTERN (GRD):
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  3
 CONVERGENCE INTERVAL (CINTERVAL):          10
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 BURN-IN ITERATIONS (NBURN):                200
 ITERATIONS (NITER):                        200
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       11456
 MC SAMPLES PER SUBJECT (ISAMPLE):          2
 RANDOM SAMPLING METHOD (RANMETHOD):        3US2P
 EXPECTATION ONLY (EONLY):                  0
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
 MASS/IMP./POST. MATRIX REFRESH SETTING (MASSREST):      -1

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
   3   4   5   6   7
 THETAS THAT ARE SIGMA-LIKE:
 
 
 MONITORING OF SEARCH:

 Stochastic/Burn-in Mode
 iteration         -200 SAEMOBJ=  -1334.15451340742
 iteration         -190 SAEMOBJ=  -12364.6316654417
 iteration         -180 SAEMOBJ=  -13172.5366059112
 iteration         -170 SAEMOBJ=  -13343.3877100016
 iteration         -160 SAEMOBJ=  -13496.2744303501
 iteration         -150 SAEMOBJ=  -13591.2705553958
 iteration         -140 SAEMOBJ=  -13673.4001688361
 iteration         -130 SAEMOBJ=  -13709.5481718753
 iteration         -120 SAEMOBJ=  -13729.2352981231
 iteration         -110 SAEMOBJ=  -13779.8654676521
 iteration         -100 SAEMOBJ=  -13762.1713264767
 iteration          -90 SAEMOBJ=  -13822.1026120278
 iteration          -80 SAEMOBJ=  -13817.2038153621
 iteration          -70 SAEMOBJ=  -13833.2044686721
 iteration          -60 SAEMOBJ=  -13856.0513064525
 iteration          -50 SAEMOBJ=  -13863.4045992691
 iteration          -40 SAEMOBJ=  -13829.7500792642
 iteration          -30 SAEMOBJ=  -13944.5922068513
 iteration          -20 SAEMOBJ=  -13932.2837958188
 iteration          -10 SAEMOBJ=  -14026.0325397139
 Reduced Stochastic/Accumulation Mode
 iteration            0 SAEMOBJ=  -14111.5255577641
 iteration           10 SAEMOBJ=  -14181.8382714078
 iteration           20 SAEMOBJ=  -14196.0459997958
 iteration           30 SAEMOBJ=  -14198.8609702305
 iteration           40 SAEMOBJ=  -14198.9916455681
 iteration           50 SAEMOBJ=  -14198.3399442224
 iteration           60 SAEMOBJ=  -14198.5795980700
 iteration           70 SAEMOBJ=  -14199.4886606500
 iteration           80 SAEMOBJ=  -14199.2047271429
 iteration           90 SAEMOBJ=  -14198.6168108129
 iteration          100 SAEMOBJ=  -14198.6677908842
 iteration          110 SAEMOBJ=  -14197.1693706831
 iteration          120 SAEMOBJ=  -14198.0472038464
 iteration          130 SAEMOBJ=  -14197.5061841477
 iteration          140 SAEMOBJ=  -14196.5518990895
 iteration          150 SAEMOBJ=  -14196.1130390702
 iteration          160 SAEMOBJ=  -14196.6933694492
 iteration          170 SAEMOBJ=  -14196.3792813015
 iteration          180 SAEMOBJ=  -14196.0897241339
 iteration          190 SAEMOBJ=  -14196.2025720543
 iteration          200 SAEMOBJ=  -14196.1265244115
 
 #TERM:
 STOCHASTIC PORTION WAS NOT COMPLETED
 REDUCED STOCHASTIC PORTION WAS COMPLETED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         1.9594E-06 -2.4940E-06 -1.1360E-05 -5.0580E-05 -8.7098E-06
 SE:             2.0994E-03  1.5112E-02  3.9120E-03  8.1675E-03  1.9726E-02
 N:                     125         125         125         125         125
 
 P VAL.:         9.9926E-01  9.9987E-01  9.9768E-01  9.9506E-01  9.9965E-01
 
 ETASHRINKSD(%)  3.4521E+01  8.6075E+00  2.9888E+01  4.0591E+01  1.3087E+00
 ETASHRINKVR(%)  5.7125E+01  1.6474E+01  5.0844E+01  6.4705E+01  2.6002E+00
 EBVSHRINKSD(%)  3.4485E+01  8.6228E+00  2.9877E+01  4.0599E+01  1.3272E+00
 EBVSHRINKVR(%)  5.7077E+01  1.6502E+01  5.0828E+01  6.4715E+01  2.6368E+00
 EPSSHRINKSD(%)  1.0000E-10
 EPSSHRINKVR(%)  1.0000E-10
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         3125
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    5743.36583252920     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -14196.1265244115     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -8452.76069188232     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           625
 NIND*NETA*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1148.67316650584     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -14196.1265244115     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -13047.4533579057     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 #TERE:
 Elapsed estimation  time in seconds:   403.67
 Elapsed covariance  time in seconds:     0.05
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION               ********************
 #OBJT:**************                        FINAL VALUE OF LIKELIHOOD FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -14196.127       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION               ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         2.50E-01  1.00E+00  4.93E-01  4.97E-02  5.12E+01  1.04E+00  2.00E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        1.30E-03
 
 ETA2
+       -1.27E-03  3.45E-02
 
 ETA3
+       -2.00E-03  3.41E-04  3.92E-03
 
 ETA4
+       -2.60E-03  3.97E-03  6.28E-03  2.38E-02
 
 ETA5
+       -5.07E-03  3.10E-03  9.64E-03  1.90E-02  5.03E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        9.40E-03
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.60E-02
 
 ETA2
+       -1.91E-01  1.86E-01
 
 ETA3
+       -8.86E-01  2.93E-02  6.26E-02
 
 ETA4
+       -4.69E-01  1.39E-01  6.50E-01  1.54E-01
 
 ETA5
+       -6.27E-01  7.45E-02  6.86E-01  5.49E-01  2.24E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        9.70E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION               ********************
 ********************                          STANDARD ERROR OF ESTIMATE (S)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         0.00E+00  0.00E+00  4.96E-03  1.06E-03  7.44E-01  3.79E-02  4.70E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        1.01E-03
 
 ETA2
+        2.04E-03  6.15E-03
 
 ETA3
+        1.52E-03  2.92E-03  2.57E-03
 
 ETA4
+        3.45E-03  9.16E-03  5.73E-03  1.71E-02
 
 ETA5
+        2.20E-03  5.59E-03  3.33E-03  9.79E-03  9.25E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        3.01E-04
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        1.41E-02
 
 ETA2
+        3.18E-01  1.66E-02
 
 ETA3
+        1.83E-01  2.51E-01  2.05E-02
 
 ETA4
+        4.77E-01  2.94E-01  3.14E-01  5.53E-02
 
 ETA5
+        2.49E-01  1.32E-01  2.29E-01  2.01E-01  2.06E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.55E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION               ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (S)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 TH 1
+       .........
 
 TH 2
+       ......... .........
 
 TH 3
+        0.00E+00  0.00E+00  2.46E-05
 
 TH 4
+        0.00E+00  0.00E+00 -8.39E-07  1.12E-06
 
 TH 5
+        0.00E+00  0.00E+00 -2.82E-03  1.73E-05  5.54E-01
 
 TH 6
+        0.00E+00  0.00E+00 -9.93E-05  1.00E-05  1.90E-02  1.44E-03
 
 TH 7
+        0.00E+00  0.00E+00 -1.20E-03  8.90E-05  1.54E-01  8.73E-03  2.21E-01
 
 OM11
+        0.00E+00  0.00E+00  1.23E-06  7.64E-08 -1.19E-04 -3.84E-06 -8.39E-05  1.03E-06
 
 OM12
+        0.00E+00  0.00E+00 -1.74E-06  3.69E-07  1.09E-04  2.13E-05  5.62E-05  7.02E-08  4.16E-06
 
 OM13
+        0.00E+00  0.00E+00 -2.18E-06 -2.66E-07  2.57E-04  6.93E-06  1.28E-04 -1.33E-06 -2.18E-07  2.30E-06
 
 OM14
+        0.00E+00  0.00E+00 -5.69E-06 -8.45E-08  6.59E-04  3.59E-05  2.52E-04 -2.01E-06  1.95E-06  3.70E-06  1.19E-05
 
 OM15
+        0.00E+00  0.00E+00 -4.28E-06  2.38E-07  4.47E-04  2.20E-05  2.37E-04 -1.04E-06  1.30E-06  1.40E-06  4.34E-06  4.83E-06
 
 OM22
+        0.00E+00  0.00E+00  3.92E-06 -1.20E-07 -2.32E-04 -2.22E-05 -5.03E-04  1.20E-06 -2.12E-07 -1.40E-06 -4.95E-06 -2.19E-06
          3.78E-05
 
 OM23
+        0.00E+00  0.00E+00  1.33E-06 -1.98E-07 -2.40E-05 -2.81E-05 -7.73E-05  2.04E-08 -4.83E-06 -1.21E-07 -3.04E-06 -1.51E-06
          7.25E-07  8.52E-06
 
 OM24
+        0.00E+00  0.00E+00  1.04E-05 -7.42E-07 -1.20E-03 -1.52E-04 -8.99E-04  1.78E-06 -1.05E-05 -2.33E-06 -1.54E-05 -6.68E-06
          2.06E-05  1.89E-05  8.38E-05
 
 OM25
+        0.00E+00  0.00E+00  1.76E-06 -1.47E-06 -1.89E-04 -4.75E-05 -4.54E-04  5.77E-07 -5.57E-06 -2.78E-08 -2.82E-06 -3.21E-06
          2.36E-06  7.47E-06  2.05E-05  3.13E-05
 
 OM33
+        0.00E+00  0.00E+00  2.96E-06  4.61E-07 -3.69E-04 -6.43E-06 -1.56E-04  1.64E-06  5.98E-07 -3.44E-06 -5.29E-06 -1.70E-06
          9.39E-07 -2.91E-07  1.43E-06 -1.44E-06  6.59E-06
 
 OM34
+        0.00E+00  0.00E+00  8.80E-06 -1.57E-07 -1.11E-03 -4.47E-05 -4.65E-04  2.61E-06 -1.15E-06 -5.78E-06 -1.50E-05 -5.59E-06
          4.01E-06  2.18E-06  1.30E-05  2.54E-06  1.15E-05  3.29E-05
 
 OM35
+        0.00E+00  0.00E+00  4.88E-06 -8.09E-08 -4.57E-04 -2.11E-05 -2.75E-04  9.41E-07 -1.03E-06 -1.93E-06 -5.92E-06 -5.32E-06
          2.67E-06  1.21E-06  7.32E-06 -4.94E-07  4.01E-06  1.06E-05  1.11E-05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 OM44
+        0.00E+00  0.00E+00  2.28E-05 -4.67E-06 -2.86E-03 -2.16E-04 -1.64E-03  2.50E-06 -1.15E-05 -6.26E-06 -2.77E-05 -1.65E-05
          1.22E-05  1.73E-05  7.81E-05  2.38E-05  1.40E-05  6.63E-05  2.56E-05  2.91E-04
 
 OM45
+        0.00E+00  0.00E+00  1.12E-05 -1.86E-06 -1.17E-03 -1.07E-04 -9.55E-04  1.26E-06 -6.13E-06 -2.79E-06 -1.94E-05 -1.32E-05
          9.83E-06  8.76E-06  4.02E-05  1.14E-05  5.06E-06  2.71E-05  2.33E-05  1.08E-04  9.59E-05
 
 OM55
+        0.00E+00  0.00E+00  8.16E-06 -2.06E-06 -1.09E-03 -7.89E-05 -7.48E-04  1.04E-06 -2.02E-06 -1.99E-06 -8.76E-06 -1.06E-05
          8.02E-06  2.12E-06  1.02E-05  1.33E-05  2.44E-06  1.43E-05  1.15E-05  5.57E-05  5.17E-05  8.55E-05
 
 SG11
+        0.00E+00  0.00E+00  1.74E-07  8.63E-09 -2.79E-05 -1.59E-06 -3.89E-06 -5.83E-08 -1.38E-07  4.71E-08 -1.42E-07 -3.21E-08
          7.86E-08  1.92E-07  3.41E-07  1.02E-07 -9.31E-08  6.85E-08  7.95E-08  8.84E-08  4.58E-07 -1.06E-07  9.05E-08
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION               ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (S)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 TH 1
+       .........
 
 TH 2
+       ......... .........
 
 TH 3
+        0.00E+00  0.00E+00  4.96E-03
 
 TH 4
+        0.00E+00  0.00E+00 -1.60E-01  1.06E-03
 
 TH 5
+        0.00E+00  0.00E+00 -7.63E-01  2.19E-02  7.44E-01
 
 TH 6
+        0.00E+00  0.00E+00 -5.28E-01  2.49E-01  6.73E-01  3.79E-02
 
 TH 7
+        0.00E+00  0.00E+00 -5.15E-01  1.79E-01  4.39E-01  4.90E-01  4.70E-01
 
 OM11
+        0.00E+00  0.00E+00  2.45E-01  7.12E-02 -1.58E-01 -9.99E-02 -1.76E-01  1.01E-03
 
 OM12
+        0.00E+00  0.00E+00 -1.72E-01  1.71E-01  7.21E-02  2.75E-01  5.87E-02  3.40E-02  2.04E-03
 
 OM13
+        0.00E+00  0.00E+00 -2.90E-01 -1.66E-01  2.27E-01  1.20E-01  1.79E-01 -8.69E-01 -7.06E-02  1.52E-03
 
 OM14
+        0.00E+00  0.00E+00 -3.33E-01 -2.31E-02  2.57E-01  2.74E-01  1.56E-01 -5.75E-01  2.77E-01  7.07E-01  3.45E-03
 
 OM15
+        0.00E+00  0.00E+00 -3.93E-01  1.02E-01  2.73E-01  2.64E-01  2.29E-01 -4.68E-01  2.90E-01  4.19E-01  5.72E-01  2.20E-03
 
 OM22
+        0.00E+00  0.00E+00  1.28E-01 -1.84E-02 -5.07E-02 -9.53E-02 -1.74E-01  1.93E-01 -1.69E-02 -1.50E-01 -2.33E-01 -1.62E-01
          6.15E-03
 
 OM23
+        0.00E+00  0.00E+00  9.21E-02 -6.39E-02 -1.10E-02 -2.54E-01 -5.64E-02  6.90E-03 -8.12E-01 -2.73E-02 -3.02E-01 -2.35E-01
          4.04E-02  2.92E-03
 
 OM24
+        0.00E+00  0.00E+00  2.29E-01 -7.65E-02 -1.77E-01 -4.37E-01 -2.09E-01  1.92E-01 -5.62E-01 -1.68E-01 -4.87E-01 -3.32E-01
          3.66E-01  7.07E-01  9.16E-03
 
 OM25
+        0.00E+00  0.00E+00  6.37E-02 -2.48E-01 -4.54E-02 -2.24E-01 -1.73E-01  1.02E-01 -4.89E-01 -3.28E-03 -1.46E-01 -2.62E-01
          6.87E-02  4.58E-01  4.00E-01  5.59E-03
 
 OM33
+        0.00E+00  0.00E+00  2.32E-01  1.70E-01 -1.93E-01 -6.60E-02 -1.29E-01  6.31E-01  1.14E-01 -8.84E-01 -5.97E-01 -3.01E-01
          5.95E-02 -3.88E-02  6.07E-02 -1.00E-01  2.57E-03
 
 OM34
+        0.00E+00  0.00E+00  3.10E-01 -2.59E-02 -2.60E-01 -2.05E-01 -1.73E-01  4.49E-01 -9.85E-02 -6.65E-01 -7.57E-01 -4.43E-01
          1.14E-01  1.30E-01  2.47E-01  7.94E-02  7.82E-01  5.73E-03
 
 OM35
+        0.00E+00  0.00E+00  2.95E-01 -2.29E-02 -1.84E-01 -1.67E-01 -1.75E-01  2.79E-01 -1.51E-01 -3.82E-01 -5.14E-01 -7.26E-01
          1.30E-01  1.25E-01  2.40E-01 -2.65E-02  4.68E-01  5.53E-01  3.33E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 OM44
+        0.00E+00  0.00E+00  2.70E-01 -2.58E-01 -2.25E-01 -3.34E-01 -2.04E-01  1.45E-01 -3.31E-01 -2.42E-01 -4.70E-01 -4.40E-01
          1.16E-01  3.47E-01  5.00E-01  2.50E-01  3.20E-01  6.77E-01  4.49E-01  1.71E-02
 
 OM45
+        0.00E+00  0.00E+00  2.31E-01 -1.79E-01 -1.60E-01 -2.88E-01 -2.08E-01  1.28E-01 -3.07E-01 -1.88E-01 -5.74E-01 -6.14E-01
          1.63E-01  3.07E-01  4.48E-01  2.08E-01  2.01E-01  4.83E-01  7.15E-01  6.44E-01  9.79E-03
 
 OM55
+        0.00E+00  0.00E+00  1.78E-01 -2.10E-01 -1.59E-01 -2.25E-01 -1.72E-01  1.12E-01 -1.07E-01 -1.42E-01 -2.75E-01 -5.21E-01
          1.41E-01  7.85E-02  1.21E-01  2.56E-01  1.03E-01  2.70E-01  3.73E-01  3.53E-01  5.71E-01  9.25E-03
 
 SG11
+        0.00E+00  0.00E+00  1.17E-01  2.71E-02 -1.24E-01 -1.39E-01 -2.75E-02 -1.91E-01 -2.24E-01  1.03E-01 -1.37E-01 -4.85E-02
          4.25E-02  2.19E-01  1.24E-01  6.09E-02 -1.21E-01  3.97E-02  7.92E-02  1.72E-02  1.55E-01 -3.82E-02  3.01E-04
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION               ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (S)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 TH 1
+       .........
 
 TH 2
+       ......... .........
 
 TH 3
+        0.00E+00  0.00E+00  1.34E+05
 
 TH 4
+        0.00E+00  0.00E+00  1.04E+05  1.32E+06
 
 TH 5
+        0.00E+00  0.00E+00  6.18E+02  8.96E+02  7.04E+00
 
 TH 6
+        0.00E+00  0.00E+00 -2.66E+03 -1.39E+04 -6.54E+01  1.97E+03
 
 TH 7
+        0.00E+00  0.00E+00  2.97E+02 -8.19E+01  5.11E-01 -3.26E+01  7.53E+00
 
 OM11
+        0.00E+00  0.00E+00  7.93E+04  5.55E+05 -4.27E+02 -2.06E+03 -3.28E+02  1.01E+07
 
 OM12
+        0.00E+00  0.00E+00  5.88E+04 -1.09E+04  7.42E+01 -3.04E+03  1.43E+02  8.97E+05  1.05E+06
 
 OM13
+        0.00E+00  0.00E+00  1.53E+05  6.87E+05 -2.39E+02 -2.08E+03 -1.07E+03  1.12E+07  1.67E+06  1.79E+07
 
 OM14
+        0.00E+00  0.00E+00 -1.76E+04 -1.01E+05  2.60E+01 -2.24E+01  4.30E+02 -1.34E+06 -4.08E+05 -2.77E+06  1.01E+06
 
 OM15
+        0.00E+00  0.00E+00  4.59E+04  1.00E+05 -1.01E+02  1.34E+03  2.19E+01  9.55E+05  1.57E+04  3.10E+05 -2.68E+05  1.02E+06
 
 OM22
+        0.00E+00  0.00E+00 -2.96E+03  4.23E+03 -2.82E+01 -1.03E+03  5.73E+01 -4.89E+03  5.33E+03  3.87E+04 -6.59E+03 -8.33E+03
          3.80E+04
 
 OM23
+        0.00E+00  0.00E+00  1.31E+04 -2.36E+04 -2.16E+02 -1.52E+03 -3.48E+01  8.27E+05  5.79E+05  1.19E+06 -2.48E+05  5.96E+04
          3.73E+04  6.44E+05
 
 OM24
+        0.00E+00  0.00E+00 -5.36E+03 -4.85E+04 -2.63E+01  2.96E+03  1.17E+01 -2.07E+05 -6.76E+04 -2.94E+05  1.16E+05 -3.90E+04
         -2.08E+04 -1.15E+05  6.21E+04
 
 OM25
+        0.00E+00  0.00E+00  1.37E+04  4.58E+04  2.78E+01 -4.89E+02  8.33E+01  1.30E+04  6.73E+04  1.02E+05 -6.10E+04  6.04E+04
          2.45E+03  4.23E+03 -1.32E+04  6.02E+04
 
 OM33
+        0.00E+00  0.00E+00  5.04E+04  7.78E+04 -5.43E+01 -1.98E+02 -3.17E+02  3.27E+06  5.35E+05  6.61E+06 -1.03E+06 -1.78E+05
          3.05E+04  3.60E+05 -9.87E+04  5.03E+04  3.32E+06
 
 OM34
+        0.00E+00  0.00E+00 -1.48E+04 -4.56E+04  4.62E+01 -3.01E+02  1.32E+02 -6.78E+05 -2.00E+05 -1.43E+06  5.01E+05 -1.16E+05
         -1.06E+04 -1.21E+05  6.75E+04 -3.91E+04 -8.03E+05  4.35E+05
 
 OM35
+        0.00E+00  0.00E+00 -6.84E+03  2.18E+04 -9.49E+01 -2.53E+02 -2.52E+01  2.54E+05  7.98E+04 -4.93E+04 -1.95E+05  4.88E+05
         -6.73E+03  5.21E+04 -2.36E+04  5.34E+04 -2.59E+05 -9.13E+04  5.20E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 OM44
+        0.00E+00  0.00E+00  2.58E+03  3.33E+04  1.19E+01 -9.79E+01 -2.68E+01  1.22E+05  3.61E+04  1.91E+05 -8.41E+04  3.65E+04
          2.64E+03  2.17E+04 -1.65E+04  6.59E+03  8.46E+04 -6.93E+04  3.04E+04  2.00E+04
 
 OM45
+        0.00E+00  0.00E+00  2.31E+03 -1.34E+04 -1.68E+01 -3.76E+01  9.59E+01 -1.67E+05 -5.44E+04 -3.20E+05  1.44E+05 -6.68E+04
          2.82E+03 -2.71E+04  7.15E+03 -9.50E+03 -7.50E+04  6.39E+04 -1.11E+05 -1.91E+04  6.67E+04
 
 OM55
+        0.00E+00  0.00E+00  2.98E+03  1.85E+04  4.52E+00  7.29E+02 -2.86E+01  1.49E+05  6.84E+03  1.81E+05 -4.29E+04  6.34E+04
         -5.11E+03  5.73E+03  3.53E+03 -3.21E+03  5.20E+04 -2.07E+04  3.16E+04  4.05E+03 -2.08E+04  2.66E+04
 
 SG11
+        0.00E+00  0.00E+00 -2.41E+04 -1.10E+05  3.50E+02  1.00E+04 -4.18E+02  3.46E+06  3.59E+05  3.65E+06 -2.82E+05  1.25E+05
         -6.14E+04 -1.23E+05  4.03E+04 -1.04E+03  1.79E+06 -5.02E+05 -8.89E+03  1.09E+05 -1.57E+05  1.29E+05  1.51E+07
 
1
 
 
 #TBLN:      2
 #METH: Importance Sampling
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            960
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
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      3
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     3
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): pdlidr.ext
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
 ITERATIONS (NITER):                        200
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       11456
 MC SAMPLES PER SUBJECT (ISAMPLE):          300
 RANDOM SAMPLING METHOD (RANMETHOD):        3US2P
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
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   4
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR TABLE/SCATTER STEP:
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
   3   4   5   6   7
 THETAS THAT ARE SIGMA-LIKE:
 
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=  -10439.9393683531 eff.=     320. Smpl.=     300. Fit.= 0.97030
 iteration            1 OBJ=  -10443.1162298101 eff.=    2580. Smpl.=     300. Fit.= 0.92194
 iteration            2 OBJ=  -10424.7001305120 eff.=     576. Smpl.=     300. Fit.= 0.91449
 iteration            3 OBJ=  -10425.3366934434 eff.=    1713. Smpl.=     300. Fit.= 0.91417
 iteration            4 OBJ=  -10415.5151480657 eff.=     513. Smpl.=     300. Fit.= 0.91405
 iteration            5 OBJ=  -10435.9004217161 eff.=     381. Smpl.=     300. Fit.= 0.91992
 iteration            6 OBJ=  -10441.3593220472 eff.=     123. Smpl.=     300. Fit.= 0.92889
 iteration            7 OBJ=  -10441.4058418530 eff.=     117. Smpl.=     300. Fit.= 0.92882
 iteration            8 OBJ=  -10427.9204199874 eff.=     588. Smpl.=     300. Fit.= 0.91690
 iteration            9 OBJ=  -10421.3557069629 eff.=     475. Smpl.=     300. Fit.= 0.91432
 iteration           10 OBJ=  -10429.8222758459 eff.=     557. Smpl.=     300. Fit.= 0.92049
 iteration           11 OBJ=  -10437.1574461985 eff.=     255. Smpl.=     300. Fit.= 0.92423
 iteration           12 OBJ=  -10436.2212693874 eff.=     475. Smpl.=     300. Fit.= 0.92161
 iteration           13 OBJ=  -10434.3883547483 eff.=     393. Smpl.=     300. Fit.= 0.91824
 iteration           14 OBJ=  -10439.4381166971 eff.=     118. Smpl.=     300. Fit.= 0.92980
 iteration           15 OBJ=  -10423.5914747619 eff.=     561. Smpl.=     300. Fit.= 0.91026
 iteration           16 OBJ=  -10407.9995874935 eff.=     800. Smpl.=     300. Fit.= 0.90901
 iteration           17 OBJ=  -10424.9462850885 eff.=    1616. Smpl.=     300. Fit.= 0.91483
 iteration           18 OBJ=  -10424.8713583340 eff.=    5919. Smpl.=     300. Fit.= 0.91945
 iteration           19 OBJ=  -10410.5624019736 eff.=     619. Smpl.=     300. Fit.= 0.91005
 iteration           20 OBJ=  -10413.7132452705 eff.=    3602. Smpl.=     300. Fit.= 0.91647
 iteration           21 OBJ=  -10407.3627319124 eff.=     945. Smpl.=     300. Fit.= 0.91231
 iteration           22 OBJ=  -10418.2164677340 eff.=     662. Smpl.=     300. Fit.= 0.91172
 iteration           23 OBJ=  -10437.0099416430 eff.=     321. Smpl.=     300. Fit.= 0.92374
 iteration           24 OBJ=  -10442.0369941914 eff.=     121. Smpl.=     300. Fit.= 0.92800
 iteration           25 OBJ=  -10441.0368145446 eff.=     387. Smpl.=     300. Fit.= 0.92532
 iteration           26 OBJ=  -10428.6897280624 eff.=     528. Smpl.=     300. Fit.= 0.91396
 iteration           27 OBJ=  -10436.0342958167 eff.=     436. Smpl.=     300. Fit.= 0.91835
 iteration           28 OBJ=  -10439.5726884412 eff.=     273. Smpl.=     300. Fit.= 0.92567
 iteration           29 OBJ=  -10436.7480255144 eff.=     421. Smpl.=     300. Fit.= 0.92075
 iteration           30 OBJ=  -10439.0724744170 eff.=     181. Smpl.=     300. Fit.= 0.92510
 iteration           31 OBJ=  -10443.0445672840 eff.=     121. Smpl.=     300. Fit.= 0.92968
 iteration           32 OBJ=  -10441.2002557590 eff.=     397. Smpl.=     300. Fit.= 0.92556
 iteration           33 OBJ=  -10437.6979764095 eff.=     413. Smpl.=     300. Fit.= 0.92234
 iteration           34 OBJ=  -10434.5097059657 eff.=     542. Smpl.=     300. Fit.= 0.91807
 iteration           35 OBJ=  -10437.4169559376 eff.=    1101. Smpl.=     300. Fit.= 0.91884
 iteration           36 OBJ=  -10420.9209725408 eff.=     515. Smpl.=     300. Fit.= 0.91751
 iteration           37 OBJ=  -10427.5365337824 eff.=     585. Smpl.=     300. Fit.= 0.91985
 iteration           38 OBJ=  -10435.9177159115 eff.=     161. Smpl.=     300. Fit.= 0.92608
 iteration           39 OBJ=  -10443.1623246459 eff.=     121. Smpl.=     300. Fit.= 0.92783
 iteration           40 OBJ=  -10444.4777215410 eff.=     119. Smpl.=     300. Fit.= 0.92856
 iteration           41 OBJ=  -10444.7974124528 eff.=     119. Smpl.=     300. Fit.= 0.92870
 iteration           42 OBJ=  -10441.3136988091 eff.=     247. Smpl.=     300. Fit.= 0.92504
 iteration           43 OBJ=  -10440.5759138531 eff.=     476. Smpl.=     300. Fit.= 0.92144
 iteration           44 OBJ=  -10439.2013211037 eff.=     377. Smpl.=     300. Fit.= 0.92316
 iteration           45 OBJ=  -10435.0157792973 eff.=     485. Smpl.=     300. Fit.= 0.91778
 iteration           46 OBJ=  -10436.0865863423 eff.=     963. Smpl.=     300. Fit.= 0.91362
 iteration           47 OBJ=  -10437.1288544297 eff.=     488. Smpl.=     300. Fit.= 0.91636
 iteration           48 OBJ=  -10441.1136196937 eff.=     179. Smpl.=     300. Fit.= 0.92332
 iteration           49 OBJ=  -10442.4241123711 eff.=     233. Smpl.=     300. Fit.= 0.92601
 iteration           50 OBJ=  -10436.8678264377 eff.=     649. Smpl.=     300. Fit.= 0.91733
 iteration           51 OBJ=  -10433.1270817545 eff.=     527. Smpl.=     300. Fit.= 0.91847
 iteration           52 OBJ=  -10441.9233677004 eff.=     626. Smpl.=     300. Fit.= 0.92133
 iteration           53 OBJ=  -10437.0045455270 eff.=    1121. Smpl.=     300. Fit.= 0.91195
 iteration           54 OBJ=  -10433.0812494834 eff.=     350. Smpl.=     300. Fit.= 0.91957
 iteration           55 OBJ=  -10441.8491282960 eff.=     245. Smpl.=     300. Fit.= 0.92426
 iteration           56 OBJ=  -10438.0492235518 eff.=     315. Smpl.=     300. Fit.= 0.92077
 iteration           57 OBJ=  -10443.1642980523 eff.=     723. Smpl.=     300. Fit.= 0.92001
 iteration           58 OBJ=  -10441.8717796316 eff.=     200. Smpl.=     300. Fit.= 0.92631
 iteration           59 OBJ=  -10443.8455280523 eff.=     118. Smpl.=     300. Fit.= 0.93000
 iteration           60 OBJ=  -10427.7753834133 eff.=     343. Smpl.=     300. Fit.= 0.91983
 iteration           61 OBJ=  -10432.5983647880 eff.=     506. Smpl.=     300. Fit.= 0.91966
 iteration           62 OBJ=  -10441.1679138303 eff.=     368. Smpl.=     300. Fit.= 0.92326
 iteration           63 OBJ=  -10440.9909856948 eff.=     413. Smpl.=     300. Fit.= 0.92197
 iteration           64 OBJ=  -10443.1385222775 eff.=     163. Smpl.=     300. Fit.= 0.92595
 iteration           65 OBJ=  -10447.1507487663 eff.=     120. Smpl.=     300. Fit.= 0.92843
 iteration           66 OBJ=  -10446.0101012608 eff.=     121. Smpl.=     300. Fit.= 0.93060
 iteration           67 OBJ=  -10444.6413599127 eff.=     118. Smpl.=     300. Fit.= 0.92980
 iteration           68 OBJ=  -10426.5117208563 eff.=     470. Smpl.=     300. Fit.= 0.92000
 iteration           69 OBJ=  -10424.0627708535 eff.=     273. Smpl.=     300. Fit.= 0.91631
 iteration           70 OBJ=  -10431.5198057161 eff.=     423. Smpl.=     300. Fit.= 0.92122
 iteration           71 OBJ=  -10438.5742062781 eff.=     413. Smpl.=     300. Fit.= 0.92329
 iteration           72 OBJ=  -10444.7701903731 eff.=     130. Smpl.=     300. Fit.= 0.92920
 iteration           73 OBJ=  -10440.3969112776 eff.=     256. Smpl.=     300. Fit.= 0.92651
 iteration           74 OBJ=  -10427.4906116031 eff.=     340. Smpl.=     300. Fit.= 0.92176
 iteration           75 OBJ=  -10426.1483851696 eff.=     869. Smpl.=     300. Fit.= 0.91965
 iteration           76 OBJ=  -10419.1970078918 eff.=     501. Smpl.=     300. Fit.= 0.91970
 iteration           77 OBJ=  -10433.8945331807 eff.=     484. Smpl.=     300. Fit.= 0.91910
 iteration           78 OBJ=  -10441.9181938612 eff.=     205. Smpl.=     300. Fit.= 0.92364
 iteration           79 OBJ=  -10445.1007365949 eff.=     314. Smpl.=     300. Fit.= 0.92469
 iteration           80 OBJ=  -10444.2184101066 eff.=     173. Smpl.=     300. Fit.= 0.92640
 iteration           81 OBJ=  -10443.9719376515 eff.=     282. Smpl.=     300. Fit.= 0.92470
 iteration           82 OBJ=  -10442.7320826500 eff.=     179. Smpl.=     300. Fit.= 0.92604
 iteration           83 OBJ=  -10443.8392751119 eff.=     297. Smpl.=     300. Fit.= 0.92588
 iteration           84 OBJ=  -10440.7712465726 eff.=     282. Smpl.=     300. Fit.= 0.92142
 iteration           85 OBJ=  -10439.2064233294 eff.=     205. Smpl.=     300. Fit.= 0.92678
 iteration           86 OBJ=  -10416.1181465338 eff.=     578. Smpl.=     300. Fit.= 0.91211
 iteration           87 OBJ=  -10424.1387485405 eff.=     548. Smpl.=     300. Fit.= 0.91742
 iteration           88 OBJ=  -10443.8844065545 eff.=     189. Smpl.=     300. Fit.= 0.92641
 iteration           89 OBJ=  -10442.5427477637 eff.=     118. Smpl.=     300. Fit.= 0.93007
 iteration           90 OBJ=  -10433.3232531299 eff.=     548. Smpl.=     300. Fit.= 0.91999
 iteration           91 OBJ=  -10433.9967459100 eff.=     274. Smpl.=     300. Fit.= 0.92304
 iteration           92 OBJ=  -10443.5346437130 eff.=     153. Smpl.=     300. Fit.= 0.92856
 iteration           93 OBJ=  -10434.8809762477 eff.=     748. Smpl.=     300. Fit.= 0.92085
 iteration           94 OBJ=  -10430.8901754969 eff.=     410. Smpl.=     300. Fit.= 0.92313
 iteration           95 OBJ=  -10435.7422710435 eff.=     502. Smpl.=     300. Fit.= 0.92266
 iteration           96 OBJ=  -10434.7032194177 eff.=     616. Smpl.=     300. Fit.= 0.91963
 iteration           97 OBJ=  -10440.4927705390 eff.=     409. Smpl.=     300. Fit.= 0.92093
 iteration           98 OBJ=  -10440.3945469175 eff.=     153. Smpl.=     300. Fit.= 0.92762
 iteration           99 OBJ=  -10429.8548688602 eff.=     759. Smpl.=     300. Fit.= 0.91614
 iteration          100 OBJ=  -10434.0044534378 eff.=     711. Smpl.=     300. Fit.= 0.91746
 iteration          101 OBJ=  -10441.6365112641 eff.=     134. Smpl.=     300. Fit.= 0.93082
 iteration          102 OBJ=  -10434.0511138434 eff.=     951. Smpl.=     300. Fit.= 0.91828
 iteration          103 OBJ=  -10433.8457421919 eff.=     622. Smpl.=     300. Fit.= 0.91927
 iteration          104 OBJ=  -10444.6929473663 eff.=     135. Smpl.=     300. Fit.= 0.92823
 iteration          105 OBJ=  -10445.2259759429 eff.=     328. Smpl.=     300. Fit.= 0.92515
 iteration          106 OBJ=  -10441.0642722941 eff.=     313. Smpl.=     300. Fit.= 0.92159
 iteration          107 OBJ=  -10444.4254360349 eff.=     194. Smpl.=     300. Fit.= 0.92643
 iteration          108 OBJ=  -10442.2689135172 eff.=     202. Smpl.=     300. Fit.= 0.92566
 iteration          109 OBJ=  -10442.3134424982 eff.=     118. Smpl.=     300. Fit.= 0.93029
 iteration          110 OBJ=  -10423.9508963545 eff.=     508. Smpl.=     300. Fit.= 0.92071
 iteration          111 OBJ=  -10427.7403998348 eff.=     606. Smpl.=     300. Fit.= 0.92036
 iteration          112 OBJ=  -10440.3784925500 eff.=     185. Smpl.=     300. Fit.= 0.92757
 iteration          113 OBJ=  -10435.5356298737 eff.=     231. Smpl.=     300. Fit.= 0.92309
 iteration          114 OBJ=  -10439.2495580793 eff.=     417. Smpl.=     300. Fit.= 0.92299
 iteration          115 OBJ=  -10445.4805122067 eff.=     143. Smpl.=     300. Fit.= 0.92853
 iteration          116 OBJ=  -10442.4138208064 eff.=     119. Smpl.=     300. Fit.= 0.93101
 iteration          117 OBJ=  -10432.7875151388 eff.=     451. Smpl.=     300. Fit.= 0.91986
 iteration          118 OBJ=  -10436.9650875228 eff.=     648. Smpl.=     300. Fit.= 0.92179
 iteration          119 OBJ=  -10447.1633088240 eff.=     138. Smpl.=     300. Fit.= 0.92862
 iteration          120 OBJ=  -10447.1350074295 eff.=     118. Smpl.=     300. Fit.= 0.92864
 iteration          121 OBJ=  -10447.1626737700 eff.=     121. Smpl.=     300. Fit.= 0.92955
 iteration          122 OBJ=  -10447.3217728406 eff.=     121. Smpl.=     300. Fit.= 0.92960
 iteration          123 OBJ=  -10447.4560882765 eff.=     118. Smpl.=     300. Fit.= 0.93004
 iteration          124 OBJ=  -10443.0830932933 eff.=     217. Smpl.=     300. Fit.= 0.92647
 iteration          125 OBJ=  -10440.2867990818 eff.=     295. Smpl.=     300. Fit.= 0.92111
 iteration          126 OBJ=  -10444.4095965046 eff.=     181. Smpl.=     300. Fit.= 0.92539
 iteration          127 OBJ=  -10447.1426375751 eff.=     121. Smpl.=     300. Fit.= 0.92995
 iteration          128 OBJ=  -10443.3424923380 eff.=     228. Smpl.=     300. Fit.= 0.92461
 iteration          129 OBJ=  -10444.9797791170 eff.=     239. Smpl.=     300. Fit.= 0.92578
 iteration          130 OBJ=  -10446.8020159882 eff.=     121. Smpl.=     300. Fit.= 0.92934
 iteration          131 OBJ=  -10443.3857325462 eff.=     205. Smpl.=     300. Fit.= 0.92573
 iteration          132 OBJ=  -10440.6765096401 eff.=     391. Smpl.=     300. Fit.= 0.92124
 iteration          133 OBJ=  -10444.9293527919 eff.=     260. Smpl.=     300. Fit.= 0.92488
 iteration          134 OBJ=  -10446.8876343781 eff.=     122. Smpl.=     300. Fit.= 0.93036
 iteration          135 OBJ=  -10443.1510977831 eff.=     253. Smpl.=     300. Fit.= 0.92500
 iteration          136 OBJ=  -10442.2114522158 eff.=     204. Smpl.=     300. Fit.= 0.92680
 iteration          137 OBJ=  -10447.5255843739 eff.=    1184. Smpl.=     300. Fit.= 0.92501
 iteration          138 OBJ=  -10438.1052165873 eff.=     179. Smpl.=     300. Fit.= 0.92633
 iteration          139 OBJ=  -10446.9525441544 eff.=     126. Smpl.=     300. Fit.= 0.92914
 iteration          140 OBJ=  -10448.7097233043 eff.=     119. Smpl.=     300. Fit.= 0.93049
 iteration          141 OBJ=  -10420.7788492271 eff.=     689. Smpl.=     300. Fit.= 0.92194
 iteration          142 OBJ=  -10420.6087251771 eff.=     268. Smpl.=     300. Fit.= 0.92374
 iteration          143 OBJ=  -10436.3317741995 eff.=     187. Smpl.=     300. Fit.= 0.92637
 iteration          144 OBJ=  -10446.3272950137 eff.=     129. Smpl.=     300. Fit.= 0.92866
 iteration          145 OBJ=  -10446.8595727438 eff.=     120. Smpl.=     300. Fit.= 0.92846
 iteration          146 OBJ=  -10447.1679707404 eff.=     118. Smpl.=     300. Fit.= 0.92961
 iteration          147 OBJ=  -10445.2390673266 eff.=     366. Smpl.=     300. Fit.= 0.92508
 iteration          148 OBJ=  -10443.4076810136 eff.=     158. Smpl.=     300. Fit.= 0.92651
 iteration          149 OBJ=  -10447.0747041424 eff.=     121. Smpl.=     300. Fit.= 0.92855
 iteration          150 OBJ=  -10447.3398123313 eff.=     119. Smpl.=     300. Fit.= 0.92956
 iteration          151 OBJ=  -10446.6316168093 eff.=     694. Smpl.=     300. Fit.= 0.92507
 iteration          152 OBJ=  -10443.5991485499 eff.=     118. Smpl.=     300. Fit.= 0.92970
 iteration          153 OBJ=  -10436.6505147624 eff.=     278. Smpl.=     300. Fit.= 0.92379
 iteration          154 OBJ=  -10440.5712632926 eff.=     636. Smpl.=     300. Fit.= 0.92430
 iteration          155 OBJ=  -10444.4124741634 eff.=     144. Smpl.=     300. Fit.= 0.92878
 iteration          156 OBJ=  -10447.5629114631 eff.=     119. Smpl.=     300. Fit.= 0.92888
 iteration          157 OBJ=  -10447.6610143701 eff.=     121. Smpl.=     300. Fit.= 0.92920
 iteration          158 OBJ=  -10447.5036565291 eff.=     120. Smpl.=     300. Fit.= 0.92913
 iteration          159 OBJ=  -10447.2279080333 eff.=     119. Smpl.=     300. Fit.= 0.92886
 iteration          160 OBJ=  -10446.9922031049 eff.=     122. Smpl.=     300. Fit.= 0.92980
 Convergence achieved
 iteration          160 OBJ=  -10447.4057157104 eff.=     119. Smpl.=     300. Fit.= 0.92828
 
 #TERM:
 OPTIMIZATION WAS COMPLETED


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -1.8644E-06 -1.0617E-04 -3.1335E-05 -1.8045E-04  1.6924E-06
 SE:             3.1077E-03  1.4756E-02  5.6419E-03  1.0823E-02  1.9737E-02
 N:                     125         125         125         125         125
 
 P VAL.:         9.9952E-01  9.9426E-01  9.9557E-01  9.8670E-01  9.9993E-01
 
 ETASHRINKSD(%)  4.0470E+01  9.3094E+00  4.1320E+01  4.3417E+01  2.0628E+00
 ETASHRINKVR(%)  6.4562E+01  1.7752E+01  6.5567E+01  6.7983E+01  4.0830E+00
 EBVSHRINKSD(%)  4.0666E+01  9.3094E+00  4.1569E+01  4.3899E+01  2.0510E+00
 EBVSHRINKVR(%)  6.4795E+01  1.7752E+01  6.5858E+01  6.8527E+01  4.0598E+00
 EPSSHRINKSD(%)  4.6525E+00
 EPSSHRINKVR(%)  9.0885E+00
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         3125
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    5743.36583252920     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -10447.4057157104     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -4704.03988318124     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           625
  
 #TERE:
 Elapsed estimation  time in seconds:  1879.87

 Number of Negative Eigenvalues in Matrix=           2
 Most negative value=  -2597425.40671364
 Most positive value=   108589158.339712
 Forcing positive definiteness
 Root mean square deviation of matrix from original=   4.781773171969263E-002

 Elapsed covariance  time in seconds:    12.59
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -10447.406       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         2.50E-01  1.00E+00  4.91E-01  4.98E-02  5.12E+01  1.03E+00  2.00E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.43E-03
 
 ETA2
+       -2.34E-04  3.34E-02
 
 ETA3
+       -6.16E-03 -3.04E-04  1.16E-02
 
 ETA4
+       -9.57E-03  3.28E-03  2.00E-02  4.61E-02
 
 ETA5
+       -5.80E-03  2.58E-03  9.52E-03  1.93E-02  5.12E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        9.30E-03
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        5.86E-02
 
 ETA2
+       -2.19E-02  1.83E-01
 
 ETA3
+       -9.74E-01 -1.54E-02  1.08E-01
 
 ETA4
+       -7.60E-01  8.35E-02  8.63E-01  2.15E-01
 
 ETA5
+       -4.38E-01  6.24E-02  3.90E-01  3.98E-01  2.26E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        9.65E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                          STANDARD ERROR OF ESTIMATE (R)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         0.00E+00  0.00E+00  1.95E-02  1.93E-03  9.31E-01  1.20E-01  9.79E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        7.31E-03
 
 ETA2
+        8.50E-03  9.36E-03
 
 ETA3
+        8.90E-03  5.78E-03  7.70E-03
 
 ETA4
+        2.16E-02  3.17E-02  2.49E-02  6.37E-02
 
 ETA5
+        1.46E-02  1.38E-02  1.08E-02  3.50E-02  2.47E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        2.60E-04
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        6.24E-02
 
 ETA2
+        7.70E-01  2.56E-02
 
 ETA3
+        1.01E-01  2.98E-01  3.57E-02
 
 ETA4
+        4.24E-01  7.45E-01  2.37E-01  1.48E-01
 
 ETA5
+        5.54E-01  3.13E-01  2.68E-01  3.73E-01  5.45E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.35E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (R)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 TH 1
+       .........
 
 TH 2
+       ......... .........
 
 TH 3
+        0.00E+00  0.00E+00  3.82E-04
 
 TH 4
+        0.00E+00  0.00E+00 -3.24E-05  3.73E-06
 
 TH 5
+        0.00E+00  0.00E+00 -1.02E-03 -2.46E-04  8.66E-01
 
 TH 6
+        0.00E+00  0.00E+00 -2.27E-03  2.00E-04  9.48E-03  1.45E-02
 
 TH 7
+        0.00E+00  0.00E+00 -1.79E-02  1.58E-03 -1.37E-02  1.09E-01  9.58E-01
 
 OM11
+        0.00E+00  0.00E+00  1.35E-04 -1.21E-05  1.04E-03 -8.19E-04 -6.46E-03  5.35E-05
 
 OM12
+        0.00E+00  0.00E+00 -1.54E-04  1.39E-05 -1.09E-03  9.36E-04  7.38E-03 -5.88E-05  7.23E-05
 
 OM13
+        0.00E+00  0.00E+00 -1.62E-04  1.45E-05 -1.23E-03  9.82E-04  7.73E-03 -6.45E-05  7.04E-05  7.92E-05
 
 OM14
+        0.00E+00  0.00E+00 -4.00E-04  3.60E-05 -3.00E-03  2.42E-03  1.91E-02 -1.55E-04  1.75E-04  1.89E-04  4.65E-04
 
 OM15
+        0.00E+00  0.00E+00 -2.73E-04  2.44E-05 -2.05E-03  1.65E-03  1.30E-02 -1.05E-04  1.19E-04  1.26E-04  3.09E-04  2.14E-04
 
 OM22
+        0.00E+00  0.00E+00  1.46E-04 -1.32E-05  1.11E-03 -8.87E-04 -6.98E-03  5.53E-05 -6.52E-05 -6.62E-05 -1.64E-04 -1.12E-04
          8.76E-05
 
 OM23
+        0.00E+00  0.00E+00  8.96E-05 -7.96E-06  4.34E-04 -5.41E-04 -4.25E-03  3.37E-05 -4.43E-05 -4.01E-05 -1.01E-04 -6.80E-05
          3.59E-05  3.35E-05
 
 OM24
+        0.00E+00  0.00E+00  5.79E-04 -5.23E-05  4.04E-03 -3.53E-03 -2.77E-02  2.20E-04 -2.64E-04 -2.64E-04 -6.56E-04 -4.44E-04
          2.46E-04  1.64E-04  1.01E-03
 
 OM25
+        0.00E+00  0.00E+00  2.45E-04 -2.20E-05  1.90E-03 -1.48E-03 -1.17E-02  9.33E-05 -1.12E-04 -1.12E-04 -2.77E-04 -1.89E-04
          1.06E-04  6.60E-05  4.16E-04  1.91E-04
 
 OM33
+        0.00E+00  0.00E+00  1.13E-04 -1.03E-05  7.31E-04 -7.02E-04 -5.46E-03  4.83E-05 -5.02E-05 -6.26E-05 -1.45E-04 -9.19E-05
          4.70E-05  2.78E-05  1.89E-04  8.00E-05  5.93E-05
 
 OM34
+        0.00E+00  0.00E+00  4.40E-04 -3.97E-05  3.13E-03 -2.69E-03 -2.11E-02  1.74E-04 -1.95E-04 -2.16E-04 -5.29E-04 -3.45E-04
          1.81E-04  1.12E-04  7.32E-04  3.09E-04  1.77E-04  6.20E-04
 
 OM35
+        0.00E+00  0.00E+00  1.89E-04 -1.68E-05  1.19E-03 -1.15E-03 -9.00E-03  7.37E-05 -8.22E-05 -8.96E-05 -2.20E-04 -1.53E-04
          7.65E-05  4.83E-05  3.10E-04  1.30E-04  6.77E-05  2.48E-04  1.18E-04
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 OM44
+        0.00E+00  0.00E+00  1.15E-03 -1.03E-04  8.41E-03 -6.98E-03 -5.49E-02  4.44E-04 -5.09E-04 -5.43E-04 -1.35E-03 -8.91E-04
          4.75E-04  2.96E-04  1.92E-03  8.08E-04  4.21E-04  1.56E-03  6.37E-04  4.06E-03
 
 OM45
+        0.00E+00  0.00E+00  6.42E-04 -5.72E-05  4.18E-03 -3.89E-03 -3.06E-02  2.45E-04 -2.80E-04 -2.96E-04 -7.33E-04 -5.01E-04
          2.62E-04  1.63E-04  1.05E-03  4.48E-04  2.18E-04  8.26E-04  3.68E-04  2.15E-03  1.23E-03
 
 OM55
+        0.00E+00  0.00E+00  4.46E-04 -3.99E-05  3.43E-03 -2.70E-03 -2.13E-02  1.70E-04 -1.94E-04 -2.04E-04 -5.03E-04 -3.50E-04
          1.83E-04  1.11E-04  7.26E-04  3.13E-04  1.48E-04  5.61E-04  2.48E-04  1.46E-03  8.29E-04  6.09E-04
 
 SG11
+        0.00E+00  0.00E+00 -1.10E-06  1.01E-07 -9.02E-06  6.85E-06  5.34E-05 -4.72E-07  4.92E-07  5.85E-07  1.34E-06  8.82E-07
         -5.03E-07 -2.70E-07 -1.89E-06 -8.07E-07 -4.94E-07 -1.56E-06 -5.97E-07 -4.04E-06 -2.06E-06 -1.45E-06  6.77E-08
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (R)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 TH 1
+       .........
 
 TH 2
+       ......... .........
 
 TH 3
+        0.00E+00  0.00E+00  1.95E-02
 
 TH 4
+        0.00E+00  0.00E+00 -8.59E-01  1.93E-03
 
 TH 5
+        0.00E+00  0.00E+00 -5.62E-02 -1.37E-01  9.31E-01
 
 TH 6
+        0.00E+00  0.00E+00 -9.64E-01  8.60E-01  8.46E-02  1.20E-01
 
 TH 7
+        0.00E+00  0.00E+00 -9.33E-01  8.34E-01 -1.50E-02  9.22E-01  9.79E-01
 
 OM11
+        0.00E+00  0.00E+00  9.46E-01 -8.60E-01  1.53E-01 -9.31E-01 -9.03E-01  7.31E-03
 
 OM12
+        0.00E+00  0.00E+00 -9.30E-01  8.47E-01 -1.37E-01  9.15E-01  8.87E-01 -9.45E-01  8.50E-03
 
 OM13
+        0.00E+00  0.00E+00 -9.29E-01  8.47E-01 -1.48E-01  9.17E-01  8.88E-01 -9.92E-01  9.30E-01  8.90E-03
 
 OM14
+        0.00E+00  0.00E+00 -9.49E-01  8.64E-01 -1.49E-01  9.33E-01  9.05E-01 -9.84E-01  9.54E-01  9.86E-01  2.16E-02
 
 OM15
+        0.00E+00  0.00E+00 -9.54E-01  8.65E-01 -1.51E-01  9.38E-01  9.10E-01 -9.79E-01  9.53E-01  9.69E-01  9.81E-01  1.46E-02
 
 OM22
+        0.00E+00  0.00E+00  7.96E-01 -7.33E-01  1.28E-01 -7.88E-01 -7.62E-01  8.08E-01 -8.19E-01 -7.95E-01 -8.12E-01 -8.15E-01
          9.36E-03
 
 OM23
+        0.00E+00  0.00E+00  7.93E-01 -7.13E-01  8.06E-02 -7.78E-01 -7.50E-01  7.97E-01 -9.01E-01 -7.79E-01 -8.08E-01 -8.04E-01
          6.63E-01  5.78E-03
 
 OM24
+        0.00E+00  0.00E+00  9.34E-01 -8.53E-01  1.37E-01 -9.24E-01 -8.92E-01  9.49E-01 -9.77E-01 -9.35E-01 -9.58E-01 -9.56E-01
          8.28E-01  8.91E-01  3.17E-02
 
 OM25
+        0.00E+00  0.00E+00  9.07E-01 -8.24E-01  1.48E-01 -8.93E-01 -8.67E-01  9.24E-01 -9.52E-01 -9.09E-01 -9.30E-01 -9.35E-01
          8.16E-01  8.26E-01  9.50E-01  1.38E-02
 
 OM33
+        0.00E+00  0.00E+00  7.53E-01 -6.90E-01  1.02E-01 -7.58E-01 -7.25E-01  8.57E-01 -7.67E-01 -9.13E-01 -8.72E-01 -8.15E-01
          6.52E-01  6.24E-01  7.72E-01  7.53E-01  7.70E-03
 
 OM34
+        0.00E+00  0.00E+00  9.04E-01 -8.26E-01  1.35E-01 -8.97E-01 -8.65E-01  9.55E-01 -9.20E-01 -9.73E-01 -9.85E-01 -9.47E-01
          7.78E-01  7.76E-01  9.26E-01  8.97E-01  9.21E-01  2.49E-02
 
 OM35
+        0.00E+00  0.00E+00  8.92E-01 -8.03E-01  1.18E-01 -8.82E-01 -8.48E-01  9.29E-01 -8.91E-01 -9.28E-01 -9.39E-01 -9.61E-01
          7.53E-01  7.69E-01  9.01E-01  8.68E-01  8.10E-01  9.19E-01  1.08E-02
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 OM44
+        0.00E+00  0.00E+00  9.22E-01 -8.41E-01  1.42E-01 -9.10E-01 -8.81E-01  9.54E-01 -9.40E-01 -9.58E-01 -9.82E-01 -9.56E-01
          7.96E-01  8.04E-01  9.52E-01  9.18E-01  8.57E-01  9.85E-01  9.21E-01  6.37E-02
 
 OM45
+        0.00E+00  0.00E+00  9.37E-01 -8.46E-01  1.28E-01 -9.24E-01 -8.92E-01  9.57E-01 -9.39E-01 -9.49E-01 -9.71E-01 -9.79E-01
          7.99E-01  8.03E-01  9.49E-01  9.27E-01  8.08E-01  9.47E-01  9.68E-01  9.63E-01  3.50E-02
 
 OM55
+        0.00E+00  0.00E+00  9.24E-01 -8.38E-01  1.50E-01 -9.08E-01 -8.82E-01  9.43E-01 -9.24E-01 -9.30E-01 -9.45E-01 -9.68E-01
          7.93E-01  7.75E-01  9.27E-01  9.19E-01  7.77E-01  9.13E-01  9.25E-01  9.26E-01  9.59E-01  2.47E-02
 
 SG11
+        0.00E+00  0.00E+00 -2.17E-01  2.02E-01 -3.73E-02  2.19E-01  2.10E-01 -2.48E-01  2.22E-01  2.53E-01  2.39E-01  2.32E-01
         -2.07E-01 -1.80E-01 -2.29E-01 -2.25E-01 -2.47E-01 -2.41E-01 -2.11E-01 -2.44E-01 -2.26E-01 -2.26E-01  2.60E-04
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                               IMPORTANCE SAMPLING                              ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (R)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 TH 1
+       .........
 
 TH 2
+       ......... .........
 
 TH 3
+        0.00E+00  0.00E+00  1.09E+05
 
 TH 4
+        0.00E+00  0.00E+00  5.90E+04  1.31E+06
 
 TH 5
+        0.00E+00  0.00E+00  5.17E+02  6.57E+02  5.46E+00
 
 TH 6
+        0.00E+00  0.00E+00 -1.61E+03 -1.21E+04 -6.11E+01  1.79E+03
 
 TH 7
+        0.00E+00  0.00E+00  2.97E+02 -1.26E+02  1.15E+00 -2.80E+01  7.09E+00
 
 OM11
+        0.00E+00  0.00E+00  8.03E+04  3.68E+04  2.76E+03 -4.94E+04 -2.23E+02  2.39E+07
 
 OM12
+        0.00E+00  0.00E+00  3.13E+04 -2.49E+04  2.18E+02 -1.48E+03 -8.72E-01  1.09E+06  5.78E+05
 
 OM13
+        0.00E+00  0.00E+00  6.74E+05  1.47E+05  7.60E+03 -9.80E+04 -7.30E+00  4.08E+07  1.28E+06  6.81E+07
 
 OM14
+        0.00E+00  0.00E+00 -2.34E+05 -4.39E+04 -2.02E+03  2.07E+04 -5.12E+01 -9.85E+06 -1.18E+05 -1.47E+07  2.67E+06
 
 OM15
+        0.00E+00  0.00E+00  1.18E+04  2.33E+04  1.89E+02 -4.43E+03 -3.09E+01 -4.53E+05 -6.75E+04 -2.34E+05 -3.95E+05  4.81E+05
 
 OM22
+        0.00E+00  0.00E+00 -7.44E+02  3.05E+03 -4.00E-01  4.66E+00 -4.06E+00 -1.22E+04  2.54E+04 -2.01E+04  5.16E+03 -3.22E+03
          4.19E+04
 
 OM23
+        0.00E+00  0.00E+00  2.55E+04 -9.70E+03  2.67E+02 -3.51E+03  5.04E+00  4.53E+05  2.27E+05  5.01E+05  4.93E+03 -6.04E+04
          3.31E+04  2.93E+05
 
 OM24
+        0.00E+00  0.00E+00 -3.58E+03  2.42E+03 -6.76E+01  1.73E+03 -4.52E+00 -1.44E+04  2.43E+04  1.31E+04 -8.73E+03  4.65E+03
         -1.04E+04 -5.50E+04  4.48E+04
 
 OM25
+        0.00E+00  0.00E+00  2.55E+03 -9.77E+03  1.24E+01 -4.38E+02  3.81E+00 -3.49E+04  6.34E+04 -4.18E+04 -7.79E+03  2.10E+04
         -3.61E+03  2.05E+04 -1.26E+04  7.45E+04
 
 OM33
+        0.00E+00  0.00E+00  4.38E+05  8.50E+04  4.02E+03 -4.52E+04  1.12E+02  1.68E+07  3.22E+05  2.73E+07 -5.38E+06  1.12E+05
         -6.33E+03  1.48E+05  1.11E+01 -1.77E+04  1.10E+07
 
 OM34
+        0.00E+00  0.00E+00 -2.29E+05 -4.22E+04 -1.87E+03  1.91E+04 -7.55E+01 -7.28E+06 -3.83E+02 -1.08E+07  1.80E+06 -2.84E+05
          2.28E+03  2.35E+03  1.19E+04 -8.43E+02 -4.16E+06  1.46E+06
 
 OM35
+        0.00E+00  0.00E+00  2.19E+04  1.27E+04  2.14E+02 -3.60E+03 -2.44E+01 -7.74E+04 -2.80E+04  4.27E+05 -3.53E+05  1.83E+05
         -2.96E+03 -4.25E+04  5.09E+03  3.02E+04  3.21E+05 -2.89E+05  2.42E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM14      OM15  
             OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44      OM45      OM55      SG11  
 
 OM44
+        0.00E+00  0.00E+00  2.47E+04  6.10E+03  1.95E+02 -2.18E+03  1.07E+01  6.32E+05 -7.63E+03  8.93E+05 -1.19E+05  2.76E+04
          1.09E+03  5.58E+03 -9.48E+03  1.49E+03  3.59E+05 -1.42E+05  3.69E+04  2.73E+04
 
 OM45
+        0.00E+00  0.00E+00 -1.04E+04 -8.38E+03 -6.99E+01  1.20E+03 -6.36E+00 -1.69E+05 -8.59E+03 -3.42E+05  1.06E+05 -1.97E+04
          1.15E+03  5.04E+03 -2.35E+03 -1.20E+04 -1.41E+05  7.80E+04 -7.48E+04 -1.43E+04  4.94E+04
 
 OM55
+        0.00E+00  0.00E+00  1.77E+02  1.82E+03 -3.23E+00 -1.15E+02  2.94E+00  4.44E+03 -4.35E+03 -5.62E+03 -1.56E+04  6.41E+04
         -4.46E+01 -2.33E+03  1.22E+03 -3.31E+03 -3.04E+03 -7.18E+03  2.64E+04  1.81E+03 -1.38E+04  3.24E+04
 
 SG11
+        0.00E+00  0.00E+00 -2.27E+03 -7.82E+03  4.03E+01 -1.14E+03  2.14E+00  8.70E+05 -5.62E+03  5.64E+05  1.35E+05  1.34E+04
          2.02E+04 -6.27E+03 -7.81E+03  2.23E+03  3.12E+05 -9.58E+04 -3.24E+04  5.20E+04 -2.42E+03  6.09E+03  1.61E+07
 
 Elapsed postprocess time in seconds:     9.32
 Elapsed finaloutput time in seconds:     0.38
 #CPUT: Total CPU Time in Seconds,     2285.508
Stop Time: 
Fri 10/12/2018 
10:01 PM
