Wed 12/06/2017 
06:34 PM
$PROB TESTING RANDOM SAMPLERS
$ABBR FUNCTION vonmises(VQI,10)
$ABBR VECTOR VV(10)
$INPUT  AMT TVAL DV
$DATA  rsampler.csv

$PRED
MUQ=theta(1)
KAPPA=THETA(2)
IF(ICALL==4) THEN
VV(1)=1.0
VV(2)=MUQ
VV(3)=KAPPA
" CALL vonmises_RNG2(2,VV)
Y=VV(1)
DV=Y
ELSE
VQI(1)=DV
VQI(2)=MUQ
VQI(3)=KAPPA
WW=VONMISES(VQI)
Y=2.0*WW
ENDIF

$THETA 0.3 2.0

$SIMULATION (567811 NORMAL) (2933012 UNIFORM) (445678 NORMAL) SUBPROBLEMS=1 
$EST METHOD=0 MAXEVAL=9999 PRINT=1 -2LL NOTHETABOUNDTEST 
$COVR
$TABLE TVAL DV NOAPPEND NOPRINT FILE=ran_VONMISES2.tab
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  1) NM-TRAN INFERS THAT THE DATA ARE SINGLE-SUBJECT.
  NONMEM RUN CANNOT BE PARALLELIZED.
             
 (WARNING  13) WITH USER-WRITTEN PRED OR $PRED, NM-TRAN CANNOT APPEND THE
 MDV DATA ITEM.
             
 (WARNING  83) FUNCTIONS ARE USED IN ABBREVIATED CODE, BUT THE $SUBROUTINES
 RECORD DOES NOT INCLUDE THE "OTHER" OPTION.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:        6 DEC 2017
Days until program expires :4561
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.4.2
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 TESTING RANDOM SAMPLERS
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:    10000
 NO. OF DATA ITEMS IN DATA SET:   4
 ID DATA ITEM IS DATA ITEM NO.:   4
 DEP VARIABLE IS DATA ITEM NO.:   3
0LABELS FOR DATA ITEMS:
 AMT TVAL DV .ID.
0FORMAT FOR DATA:
 (3E5.0,1F2.0)

 TOT. NO. OF OBS RECS:    10000
 TOT. NO. OF INDIVIDUALS:    10000
0LENGTH OF THETA:   2
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
   0.3000E+00  0.2000E+01
0SIMULATION STEP OMITTED:    NO
 ORIGINAL DATA USED ON EACH NEW SIMULATION:         NO
 SEEDS RESET ON EACH NEW SUPERSET ITERATION:        YES
0SIMULATION RANDOM METHOD SELECTED (RANMETHOD): 4U
SEED   1 RESET TO INITIAL: YES
 SOURCE   1:
   SEED1:        567811   SEED2:             0   PSEUDO-NORMAL
SEED   2 RESET TO INITIAL: YES
 SOURCE   2:
   SEED1:       2933012   SEED2:             0   PSEUDO-UNIFORM
SEED   3 RESET TO INITIAL: YES
 SOURCE   3:
   SEED1:        445678   SEED2:             0   PSEUDO-NORMAL
 NUMBER OF SUBPROBLEMS:    1
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
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
04 COLUMNS APPENDED:    NO
 PRINTED:                NO
 HEADERS:               YES
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 TVAL DV
1
 PROBLEM NO.:           1      SUBPROBLEM NO.:           1

 SIMULATION STEP PERFORMED
 SOURCE  1:
    SEED1:        567811   SEED2:             0
 SOURCE  2:
    SEED1:    1166163886   SEED2:             0
 SOURCE  3:
    SEED1:        445678   SEED2:             0
 Elapsed simulation  time in seconds:     0.06
1
 
 
 #TBLN:      1
 #METH: First Order
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           SINGLE-SUBJECT
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               SLOW
 EPS-ETA INTERACTION:                     NO
 PRED F SET TO -2 LOG LIKELIHOOD:         YES
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): ran_vonmises2.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 WISHART PRIOR DF INTERPRETATION (WISHTYPE):0
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=NPRED
 RES=NRES
 WRES=NWRES
 IWRS=NIWRES
 IPRD=NIPRED
 IRS=NIRES
 
 MONITORING OF SEARCH:

 
0ITERATION NO.:    0    OBJECTIVE VALUE:   25301.5276651643        NO. OF FUNC. EVALS.:   3
 CUMULATIVE NO. OF FUNC. EVALS.:        3
 NPARAMETR:  3.0000E-01  2.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01
 GRADIENT:  -2.9377E+02 -1.8333E+02
 
0ITERATION NO.:    1    OBJECTIVE VALUE:   25301.4037529060        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       10
 NPARAMETR:  3.0245E-01  2.0102E+00
 PARAMETER:  1.0082E-01  1.0051E-01
 GRADIENT:  -8.8751E+01  4.8163E+02
 
0ITERATION NO.:    2    OBJECTIVE VALUE:   25301.3376810313        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       17
 NPARAMETR:  3.0472E-01  2.0055E+00
 PARAMETER:  1.0157E-01  1.0028E-01
 GRADIENT:   1.0218E+02  1.7741E+02
 
0ITERATION NO.:    3    OBJECTIVE VALUE:   25301.3190935621        NO. OF FUNC. EVALS.:   4
 CUMULATIVE NO. OF FUNC. EVALS.:       21
 NPARAMETR:  3.0350E-01  2.0028E+00
 PARAMETER:  1.0117E-01  1.0014E-01
 GRADIENT:  -4.5801E-02 -1.1840E+00
 
0ITERATION NO.:    4    OBJECTIVE VALUE:   25301.3190935621        NO. OF FUNC. EVALS.:   4
 CUMULATIVE NO. OF FUNC. EVALS.:       25
 NPARAMETR:  3.0350E-01  2.0028E+00
 PARAMETER:  1.0117E-01  1.0014E-01
 GRADIENT:  -1.2782E+01 -6.6812E+01
 
0ITERATION NO.:    5    OBJECTIVE VALUE:   25301.3190935621        NO. OF FUNC. EVALS.:   0
 CUMULATIVE NO. OF FUNC. EVALS.:       25
 NPARAMETR:  3.0350E-01  2.0028E+00
 PARAMETER:  1.0117E-01  1.0014E-01
 GRADIENT:  -1.2782E+01 -6.6812E+01
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:       25
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):            0
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:   0.000000000000000E+000
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:    25301.3190935621     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       25301.3190935621     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                             0
  
 #TERE:
 Elapsed estimation  time in seconds:     2.62
 Elapsed covariance  time in seconds:     1.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    25301.319       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2     
 
         3.04E-01  2.00E+00
 
1
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2     
 
         8.44E-03  2.48E-02
 
1
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2  
 
 TH 1
+        7.11E-05
 
 TH 2
+        1.27E-06  6.15E-04
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2  
 
 TH 1
+        8.44E-03
 
 TH 2
+        6.06E-03  2.48E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2  
 
 TH 1
+        1.41E+04
 
 TH 2
+       -2.90E+01  1.63E+03
 
 Elapsed finaloutput time in seconds:     0.48
 #CPUT: Total CPU Time in Seconds,        4.976
Stop Time: 
Wed 12/06/2017 
06:34 PM
