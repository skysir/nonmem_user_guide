Wed 12/06/2017 
12:32 PM
$PROB  AD3TR4
$INPUT  C SET ID JID TIME CONC AMT RATE EVID MDV CMT DV LOQ TYPE

$DATA  ad3tr4a.csv IGNORE = C
$SUBROUTINES ADVAN3 TRANS4

$PK
MU_1 = THETA(1)
MU_2 = THETA(2)
MU_3 = THETA(3)
MU_4 = THETA(4)
CL = EXP(MU_1 + ETA(1))
V1 = EXP(MU_2 + ETA(2))
Q = EXP(MU_3 + ETA(3))
V2 = EXP(MU_4 + ETA(4))
S1=V1

$ERROR
LAQ=3.0
SD = THETA(5)
IPRED = LOG(F)
DUM = (LOQ - IPRED) / SD
CUMD = PHI(DUM)+1.0E-10
DUMA = (LAQ - IPRED) / SD
CUMDA = PHI(DUMA)-1.0E-10
IF(TYPE.EQ.2) DV_LOQ=LOQ
IF(TYPE.EQ.3) DV_LAQ=LAQ
IF (TYPE .EQ. 1.OR.NPDE_MODE==1) THEN
      F_FLAG = 0
      Y = IPRED + SD * ERR(1)
ENDIF
IF (TYPE .EQ. 2.AND.NPDE_MODE==0) THEN
      F_FLAG = 1
      Y = CUMD
      MDVRES=1
ENDIF
IF (TYPE .EQ. 3.AND.NPDE_MODE==0) THEN
      F_FLAG = 1
      Y = (1.0-CUMDA)
      MDVRES=1      
ENDIF
;IF(TYPE==3.0.AND.NPDE_MODE==1) WRITE(*,*) CUMDA

$THETA 
         1.67E+00  1.60E+00  7.76E-01  2.36E+00  2.75E-01

;INITIAL values of OMEGA
$OMEGA BLOCK(4)
0.15   ;[P]
0.01  ;[F]
0.15   ;[P]
0.01  ;[F]
0.01  ;[F]
0.15   ;[P]
0.01  ;[F]
0.01  ;[F]
0.01  ;[F]
0.15   ;[P]

;Initial value of SIGMA
$SIGMA 
1.0 FIXED   ;[P]

$EST METHOD=COND INTERACTION LAPLACE MAXEVAL=9999 NSIG=3 SIGL=10 PRINT=5 NOHABORT MCETA=10
;$COV MATRIX=R PRINT=E UNCONDITIONAL
$TABLE ID TIME DV IPRED PRED RES EPRED EIPRED EIRES EIWRES CPRED CWRES ECWRES EWRES TYPE NPDE NPD CIWRES NOAPPEND ONEHEADER 
 ESAMPLE=1000 NPDTYPE=0 FILE=ad3tr4a_loq0.TAB NOPRINT SEED=16993234
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  3) THERE MAY BE AN ERROR IN THE ABBREVIATED CODE. THE FOLLOWING
 ONE OR MORE RANDOM VARIABLES ARE DEFINED WITH "IF" STATEMENTS THAT DO NOT
 PROVIDE DEFINITIONS FOR BOTH THE "THEN" AND "ELSE" CASES. IF ALL
 CONDITIONS FAIL, THE VALUES OF THESE VARIABLES WILL BE ZERO.
  
   Y

             
 (WARNING  66) WITH "LAPLACIAN" AND "INTERACTION", "NUMERICAL" AND "SLOW"
 ARE ALSO REQUIRED ON $ESTIM RECORD. NM-TRAN HAS SUPPLIED THESE OPTIONS.
  
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
 AD3TR4
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:  14
 ID DATA ITEM IS DATA ITEM NO.:   3
 DEP VARIABLE IS DATA ITEM NO.:  12
 MDV DATA ITEM IS DATA ITEM NO.: 10
0INDICES PASSED TO SUBROUTINE PRED:
   9   5   7   8   0   0  11   0   0   0   0
0LABELS FOR DATA ITEMS:
 C SET ID JID TIME CONC AMT RATE EVID MDV CMT DV LOQ TYPE
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 IPRED
0FORMAT FOR DATA:
 (2E2.0,2E4.0,E7.0,E11.0,E4.0,4E2.0,E13.0,E5.0,E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:   5
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
   0.1670E+01  0.1600E+01  0.7760E+00  0.2360E+01  0.2750E+00
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.1500E+00
                  0.1000E-01   0.1500E+00
                  0.1000E-01   0.1000E-01   0.1500E+00
                  0.1000E-01   0.1000E-01   0.1000E-01   0.1500E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    16993234
 RANMETHOD:             3U
 MC SAMPLES (ESAMPLE):    1000
 WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE
0-- TABLE   1 --
0RECORDS ONLY:    ALL
04 COLUMNS APPENDED:    NO
 PRINTED:                NO
 HEADER:                YES
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID TIME DV IPRED PRED RES EPRED EIPRED EIRES EIWRES CPRED CWRES ECWRES EWRES TYPE NPDE NPD CIWRES
1DOUBLE PRECISION PREDPP VERSION 7.4.2

 TWO COMPARTMENT MODEL (ADVAN3)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   4
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K12)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K21)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V1, Q, V2 TO K, K12, K21 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         CENTRAL      ON         NO         YES        YES        YES
    2         PERIPH.      ON         NO         YES        NO         NO
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            5           *           *           *           *
    2            *           *           *           *           *
    3            *           -           -           -           -
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
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1
 
 
 #TBLN:      1
 #METH: Laplacian Conditional Estimation with Interaction
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               SLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    YES
 NUMERICAL 2ND DERIVATIVES:               YES
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    10
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      10
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     10
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): ad3tr4a_loq0.ext
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
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 MONITORING OF SEARCH:

 
0ITERATION NO.:    0    OBJECTIVE VALUE:   3.01009013078081        NO. OF FUNC. EVALS.:  16
 CUMULATIVE NO. OF FUNC. EVALS.:       16
 NPARAMETR:  1.6700E+00  1.6000E+00  7.7600E-01  2.3600E+00  2.7500E-01  1.5000E-01  1.0000E-02  1.0000E-02  1.0000E-02  1.5000E-01
             1.0000E-02  1.0000E-02  1.5000E-01  1.0000E-02  1.5000E-01
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.3118E+01  2.5869E+02 -3.8976E+01 -3.1625E+02  1.4931E+02 -2.3621E+01  2.8432E+00 -1.3660E+00  1.5292E+01  1.0224E+01
             1.4099E+00 -1.6979E-01  1.5975E+01 -7.4212E+00  3.7730E+00
 
0ITERATION NO.:    3    OBJECTIVE VALUE:   2.75267529261327        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:       94
 NPARAMETR:  1.6691E+00  1.5836E+00  7.7721E-01  2.3897E+00  2.7337E-01  1.5003E-01  9.9998E-03  1.0001E-02  9.9949E-03  1.4999E-01
             9.9991E-03  9.9992E-03  1.4998E-01  1.0002E-02  1.5000E-01
 PARAMETER:  9.9948E-02  9.8972E-02  1.0015E-01  1.0126E-01  9.9407E-02  1.0009E-01  9.9989E-02  1.0001E-01  9.9939E-02  9.9959E-02
             9.9994E-02  1.0000E-01  9.9937E-02  1.0003E-01  9.9985E-02
 GRADIENT:  -2.3089E+01 -3.2113E+02 -1.3366E+02  2.6901E+02 -2.2589E+02 -2.4180E+01 -2.3127E+00 -4.3775E+00  1.1640E+01  4.6682E+01
            -8.0933E+00 -6.3013E+00  7.1789E+00 -8.3077E+00  2.6405E+00
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:       94
 NO. OF SIG. DIGITS IN FINAL EST.:  5.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         4.4241E-03  2.1112E-03  9.2123E-03 -7.1798E-03
 SE:             3.8954E-02  2.9308E-02  2.5232E-02  3.0527E-02
 N:                     100         100         100         100
 
 P VAL.:         9.0958E-01  9.4258E-01  7.1503E-01  8.1406E-01
 
 ETASHRINKSD(%)  1.0000E-10  2.3942E+01  3.4519E+01  2.0781E+01
 ETASHRINKVR(%)  1.0000E-10  4.2152E+01  5.7122E+01  3.7243E+01
 EBVSHRINKSD(%)  4.8267E+00  2.3261E+01  3.2475E+01  1.7732E+01
 EBVSHRINKVR(%)  9.4204E+00  4.1111E+01  5.4404E+01  3.2319E+01
 EPSSHRINKSD(%)  3.5140E+01
 EPSSHRINKVR(%)  5.7932E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          473
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    869.315852411620     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:    2.75267529261327     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       872.068527704234     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           400
  
 #TERE:
 Elapsed estimation  time in seconds:     3.59
 Elapsed postprocess time in seconds:     1.33
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                LAPLACIAN CONDITIONAL ESTIMATION WITH INTERACTION               ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************        2.753       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                LAPLACIAN CONDITIONAL ESTIMATION WITH INTERACTION               ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5     
 
         1.67E+00  1.58E+00  7.77E-01  2.39E+00  2.73E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.50E-01
 
 ETA2
+        1.00E-02  1.50E-01
 
 ETA3
+        1.00E-02  1.00E-02  1.50E-01
 
 ETA4
+        9.99E-03  1.00E-02  1.00E-02  1.50E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        3.87E-01
 
 ETA2
+        6.67E-02  3.87E-01
 
 ETA3
+        6.67E-02  6.67E-02  3.87E-01
 
 ETA4
+        6.66E-02  6.67E-02  6.67E-02  3.87E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.45
 #CPUT: Total CPU Time in Seconds,        4.914
Stop Time: 
Wed 12/06/2017 
12:32 PM
