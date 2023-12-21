Sat 04/22/2017 
09:16 AM
;Model Desc: Two Compartment model with Clearance and 
; central volume modeled with covariates age and gender
;Project Name: nm7examples
;Project ID: NO PROJECT DESCRIPTION

$PROB RUN# example2 (from sampc)
$INPUT C SET ID JID TIME DV=CONC AMT=DOSE RATE EVID MDV CMT GNDR AGE
$DATA example2.csv IGNORE=C
$SUBROUTINES ADVAN3 TRANS4

$PK
; LCLM=log transformed clearance, male
LCLM=THETA(1)
;LCLF=log transformed clearance, female.
LCLF=THETA(2)
; CLAM=CL age slope, male
CLAM=THETA(3)
; CLAF=CL age slope, female
CLAF=THETA(4)
; LV1M=log transformed V1, male
LV1M=THETA(5)
; LV1F=log transformed V1, female
LV1F=THETA(6)
; V1AM=V1 age slope, male
V1AM=THETA(7)
; V1AF=V1 age slope, female
V1AF=THETA(8)
; LAGE=log transformed age
LAGE=DLOG(AGE)

;Mean of ETA1, the inter-subject deviation of Clearance,
; is ultimately modeled as linear function of THETA(1) to THETA(4).  
; Relating thetas to Mus by linear functions is not essential for 
; ITS, IMP, or IMPMAP methods, but is very helpful for MCMC methods 
; such as SAEM and BAYES.

MU_1=(1.0-GNDR)*(LCLM+LAGE*CLAM) + GNDR*(LCLF+LAGE*CLAF)

; Mean of ETA2, the inter-subject deviation of V1, 
; is ultimately modeled as linear function of THETA(5) to THETA(8)

MU_2=(1.0-GNDR)*(LV1M+LAGE*V1AM) + GNDR*(LV1F+LAGE*V1AF)
MU_3=THETA(9)
MU_4=THETA(10)
CL=DEXP(MU_1+ETA(1))
V1=DEXP(MU_2+ETA(2))
Q=DEXP(MU_3+ETA(3))
V2=DEXP(MU_4+ETA(4))
S1=V1

$ERROR
CALLFL=0
; Option to model the residual error coefficient in THETA(11), 
; rather than in SIGMA.
SDSL=THETA(11)
W=F*SDSL
Y = F + W*EPS(1)
IPRED=F
IWRES=(DV-F)/W

;Initial THETAs
$THETA
( 0.7 ) ;[LCLM]
( 0.7 ) ;[LCLF]
( 2 )   ;[CLAM]
( 2.0);[CLAF]
( 0.7 ) ;[LV1M]
( 0.7 ) ;[LV1F]
( 2.0 )   ;[V1AM]
( 2.0 )   ;[V1AF]
( 0.7 ) ;[MU_3]
(  0.7 );[MU_4]
( 0.3 )     ;[SDSL]

;Initial OMEGAs
$OMEGA BLOCK(4)
0.5  ;[p]
0.001  ;[f]
0.5  ;[p]
0.001 ;[f]
0.001 ;[f]
0.5  ;[p]
0.001 ;[f]
0.001 ;[f]
0.001 ;[f]
0.5 ;[p]

; SIGMA is 1.0 fixed, serves as unscaled variance for EPS(1).  
; THETA(11) takes up the residual error scaling.
$SIGMA 
(1.0 FIXED)

;Prior information is important for MCMC Bayesian analysis, 
; not necessary for maximization methods
; In this example, only the OMEGAs have a prior distribution,
; the THETAS do not.
; For Bayesian methods, it is most important for at least the 
; OMEGAs to have a prior, even an uninformative one, 
; to stabilize the analysis. Only if the number of subjects
; exceeds the OMEGA dimension number by at least 100, 
; then you may get away without priors on OMEGA for BAYES analysis.
$PRIOR NWPRI
; Prior OMEGA matrix
$OMEGAP BLOCK(4) FIX VALUES(0.01,0.0)
; Degrees of freedom to OMEGA prior matrix:
$OMEGAPD 4 FIX

; The first analysis is iterative two-stage.  
; Note that the GRD specification is THETA(11) is a 
; Sigma-like parameter.  This will allow NONMEM to make
; efficient gradient evaluations for THETA(11), which is useful 
; for later IMP,IMPMAP, and SAEM methods, but has no impact on 
; ITS and BAYES methods.

$EST METHOD=ITS INTERACTION FILE=example2.ext NITER=1000 NSIG=2 
     PRINT=5 NOABORT SIGL=8 NOPRIOR=1 CTYPE=3 GRD=TS(11)

; Results of ITS serve as initial parameters for the IMP method.

$EST METHOD=IMP INTERACTION EONLY=0 MAPITER=0 NITER=100 ISAMPLE=300 
     PRINT=1 SIGL=8

; The results of IMP are used as the initial values for the SAEM method.

$EST METHOD=SAEM NBURN=3000 NITER=2000 PRINT=10 ISAMPLE=2
     CTYPE=3 CITER=10 CALPHA=0.05

; After the SAEM method, obtain good estimates of the marginal density 
; (objective function),
; along with good estimates of the standard errors.

$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=3000 
     PRINT=1 SIGL=8 SEED=123334
     CTYPE=3 CITER=10 CALPHA=0.05

; The Bayesian analysis is performed. 

$EST METHOD=BAYES INTERACTION FILE=example2.TXT NBURN=10000 
     NITER=3000 PRINT=100 NOPRIOR=0
     CTYPE=3 CITER=10 CALPHA=0.05

; Just for old-times sake, lets see what the traditional 
; FOCE method will give us.  
; And, remember to introduce a new FILE, so its results wont 
; append to our Bayesian FILE.

$EST  METHOD=COND INTERACTION MAXEVAL=9999 FILE=example2.ext NSIG=2 
  SIGL=14 PRINT=5 NOABORT NOPRIOR=1

$COV MATRIX=R UNCONDITIONAL
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

 (MU_WARNING 26) DATA ITEM(S) USED IN DEFINITION OF MU_(S) SHOULD BE CONSTANT FOR INDIV. REC.:
  GNDR AGE
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       22 APR 2017
Days until program expires :4785
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.4.0 beta 2 (nm74b2)
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 RUN# example2 (from sampc)
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     2400
 NO. OF DATA ITEMS IN DATA SET:  13
 ID DATA ITEM IS DATA ITEM NO.:   3
 DEP VARIABLE IS DATA ITEM NO.:   6
 MDV DATA ITEM IS DATA ITEM NO.: 10
0INDICES PASSED TO SUBROUTINE PRED:
   9   5   7   8   0   0  11   0   0   0   0
0LABELS FOR DATA ITEMS:
 C SET ID JID TIME CONC DOSE RATE EVID MDV CMT GNDR AGE
0FORMAT FOR DATA:
 (2E2.0,3E4.0,E11.0,E4.0,5E2.0,E6.0)

 TOT. NO. OF OBS RECS:     2000
 TOT. NO. OF INDIVIDUALS:      400
0LENGTH OF THETA:  12
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS BLOCK FORM:
  1
  1  1
  1  1  1
  1  1  1  1
  0  0  0  0  2
  0  0  0  0  2  2
  0  0  0  0  2  2  2
  0  0  0  0  2  2  2  2
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
 -0.1000E+07     0.7000E+00     0.1000E+07
 -0.1000E+07     0.7000E+00     0.1000E+07
 -0.1000E+07     0.2000E+01     0.1000E+07
 -0.1000E+07     0.2000E+01     0.1000E+07
 -0.1000E+07     0.7000E+00     0.1000E+07
 -0.1000E+07     0.7000E+00     0.1000E+07
 -0.1000E+07     0.2000E+01     0.1000E+07
 -0.1000E+07     0.2000E+01     0.1000E+07
 -0.1000E+07     0.7000E+00     0.1000E+07
 -0.1000E+07     0.7000E+00     0.1000E+07
 -0.1000E+07     0.3000E+00     0.1000E+07
  0.4000E+01     0.4000E+01     0.4000E+01
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.5000E+00
                  0.1000E-02   0.5000E+00
                  0.1000E-02   0.1000E-02   0.5000E+00
                  0.1000E-02   0.1000E-02   0.1000E-02   0.5000E+00
        2                                                                                  YES
                  0.1000E-01
                  0.0000E+00   0.1000E-01
                  0.0000E+00   0.0000E+00   0.1000E-01
                  0.0000E+00   0.0000E+00   0.0000E+00   0.1000E-01
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
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
0
 PRIOR SUBROUTINE USER-SUPPLIED
1DOUBLE PRECISION PREDPP VERSION 7.4.0 beta 2 (nm74b2)

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
0DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 OTHERWISE, ERROR SUBROUTINE CALLED ONLY WITH OBSERVATION EVENTS.
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
 NO. OF FUNCT. EVALS. ALLOWED:            2208
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      8
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     8
 NOPRIOR SETTING (NOPRIOR):                 ON
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example2.ext
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
 GRADIENT/GIBBS PATTERN (GRD):              DDDDDDDDDDS
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  3
 CONVERGENCE INTERVAL (CINTERVAL):          5
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 ITERATIONS (NITER):                        1000
 ANEAL SETTING (CONSTRAIN):                 1

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8   9  10
 THETAS THAT ARE SIGMA-LIKE:
  11
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=   43391.7048827429
 iteration            5 OBJ=  -10716.4698423578
 iteration           10 OBJ=  -10763.0072895072
 iteration           15 OBJ=  -10768.7319952295
 iteration           20 OBJ=  -10770.3950067760
 iteration           25 OBJ=  -10771.0764427879
 iteration           30 OBJ=  -10771.4099657337
 iteration           35 OBJ=  -10771.5905115024
 iteration           40 OBJ=  -10771.6945297974
 iteration           45 OBJ=  -10771.7569757847
 iteration           50 OBJ=  -10771.7955238455
 iteration           55 OBJ=  -10771.8197601345
 iteration           60 OBJ=  -10771.8351639668
 iteration           65 OBJ=  -10771.8449946358
 iteration           70 OBJ=  -10771.8512502395
 iteration           75 OBJ=  -10771.8551890581
 iteration           80 OBJ=  -10771.8576174860
 iteration           85 OBJ=  -10771.8590574267
 iteration           90 OBJ=  -10771.8598541308
 iteration           95 OBJ=  -10771.8602353296
 iteration          100 OBJ=  -10771.8603529991
 iteration          105 OBJ=  -10771.8603078827
 iteration          110 OBJ=  -10771.8601670763
 iteration          115 OBJ=  -10771.8599734420
 iteration          120 OBJ=  -10771.8597562867
 iteration          125 OBJ=  -10771.8595346484
 iteration          130 OBJ=  -10771.8593187286
 iteration          135 OBJ=  -10771.8591155198
 iteration          140 OBJ=  -10771.8589280177
 iteration          145 OBJ=  -10771.8587585453
 iteration          150 OBJ=  -10771.8586069634
 iteration          155 OBJ=  -10771.8584719967
 iteration          160 OBJ=  -10771.8583534231
 iteration          165 OBJ=  -10771.8582503632
 iteration          170 OBJ=  -10771.8581590668
 iteration          175 OBJ=  -10771.8580806407
 iteration          180 OBJ=  -10771.8580121766
 iteration          185 OBJ=  -10771.8579527504
 iteration          190 OBJ=  -10771.8579022138
 iteration          195 OBJ=  -10771.8578587179
 iteration          200 OBJ=  -10771.8578210865
 Convergence achieved
 
 #TERM:
 OPTIMIZATION WAS COMPLETED


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -4.3576E-08 -5.3097E-08 -9.8350E-08 -8.8083E-08
 SE:             4.7132E-03  2.9883E-03  2.9421E-03  3.6914E-03
 N:                     400         400         400         400
 
 P VAL.:         9.9999E-01  9.9999E-01  9.9997E-01  9.9998E-01
 
 ETASHRINKSD(%)  6.9151E+00  3.2911E+01  4.1052E+01  2.5073E+01
 ETASHRINKVR(%)  1.3352E+01  5.4991E+01  6.5251E+01  4.3860E+01
 EBVSHRINKSD(%)  6.9151E+00  3.2911E+01  4.1050E+01  2.5072E+01
 EBVSHRINKVR(%)  1.3352E+01  5.4991E+01  6.5249E+01  4.3859E+01
 EPSSHRINKSD(%)  2.6187E+01
 EPSSHRINKVR(%)  4.5517E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2000
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    3675.75413281869     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -10771.8578210865     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -7096.10368826777     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1600
  
 #TERE:
 Elapsed estimation  time in seconds:    25.18
 Elapsed covariance  time in seconds:     0.91
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -10771.858       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.30E+00  3.26E+00 -6.11E-01 -2.08E-01  7.29E-01  1.14E+00  3.37E-01  1.92E-01  6.92E-01  2.30E+00  9.99E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.03E-02
 
 ETA2
+        1.66E-04  7.96E-03
 
 ETA3
+        1.20E-03 -3.61E-04  9.99E-03
 
 ETA4
+       -6.31E-04  4.54E-04  2.02E-03  9.73E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.01E-01
 
 ETA2
+        1.83E-02  8.92E-02
 
 ETA3
+        1.18E-01 -4.04E-02  9.99E-02
 
 ETA4
+       -6.31E-02  5.16E-02  2.05E-01  9.87E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                          STANDARD ERROR OF ESTIMATE (S)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.80E-02  2.88E-02  1.11E-02  8.42E-03  4.84E-02  4.10E-02  1.34E-02  1.19E-02  1.05E-02  8.94E-03  2.98E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.01E-03
 
 ETA2
+        8.65E-04  1.50E-03
 
 ETA3
+        1.29E-03  1.44E-03  2.83E-03
 
 ETA4
+        1.12E-03  1.16E-03  1.98E-03  2.00E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        4.96E-03
 
 ETA2
+        9.47E-02  8.39E-03
 
 ETA3
+        1.20E-01  1.65E-01  1.41E-02
 
 ETA4
+        1.16E-01  1.28E-01  1.62E-01  1.01E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (S)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.45E-03
 
 TH 2
+        3.11E-06  8.32E-04
 
 TH 3
+       -4.12E-04  1.28E-06  1.23E-04
 
 TH 4
+       -4.28E-07 -2.32E-04 -1.24E-07  7.09E-05
 
 TH 5
+        6.42E-04  5.14E-05 -1.67E-04 -4.20E-06  2.34E-03
 
 TH 6
+       -1.87E-05  2.13E-04  6.23E-06 -5.21E-05  2.24E-05  1.68E-03
 
 TH 7
+       -1.68E-04 -1.53E-05  4.40E-05  1.59E-06 -6.36E-04 -9.83E-06  1.80E-04
 
 TH 8
+        6.13E-06 -5.25E-05 -1.63E-06  1.48E-05  4.11E-06 -4.71E-04 -1.00E-07  1.41E-04
 
 TH 9
+        1.75E-05  6.00E-05 -8.34E-08 -9.16E-06  1.16E-04  5.07E-05 -3.25E-05 -7.05E-06  1.11E-04
 
 TH10
+        1.38E-05  3.68E-05 -6.26E-08 -4.78E-06  1.07E-04  4.82E-05 -2.97E-05 -5.43E-06  6.60E-05  7.99E-05
 
 TH11
+       -7.59E-06 -3.65E-06  2.68E-06  1.29E-06 -8.74E-06  7.42E-07  2.49E-06 -8.00E-08  1.73E-06  1.46E-06  8.91E-06
 
 OM11
+        4.14E-06 -1.54E-07 -1.33E-06 -1.73E-07  1.28E-06  2.50E-07 -4.58E-07 -4.17E-08 -3.97E-07  1.23E-07 -4.26E-07  1.01E-06
 
 OM12
+        1.54E-06  1.63E-07 -6.67E-07 -5.19E-08  1.89E-06  7.42E-07 -3.75E-07  1.41E-08 -1.61E-07  3.35E-07 -5.62E-07  3.68E-07
          7.48E-07
 
 OM13
+        4.69E-07  1.83E-06 -4.22E-07 -5.81E-07  2.00E-06  2.41E-08 -5.61E-07  5.51E-08  7.41E-07  6.57E-07 -5.16E-07  5.33E-07
          3.32E-07  1.67E-06
 
 OM14
+        1.94E-06 -4.59E-07 -7.43E-07  2.95E-07  6.08E-06 -1.89E-07 -1.52E-06  1.38E-07  7.56E-07  9.22E-07 -3.41E-07  4.32E-07
          2.96E-07  1.03E-06  1.25E-06
 
 OM22
+        3.82E-06  6.75E-07 -1.07E-06 -9.10E-08  5.40E-06 -7.12E-06 -7.31E-07  2.13E-06 -3.58E-08 -4.58E-08 -1.07E-06  9.10E-08
          4.66E-07  6.41E-08  1.49E-07  2.24E-06
 
 OM23
+        3.10E-07  5.59E-07 -2.35E-07 -2.86E-07 -4.52E-07  1.26E-06  1.85E-07 -6.52E-08 -6.99E-07 -3.54E-07 -6.87E-07  2.99E-07
          4.27E-07  7.84E-07  5.05E-07  5.30E-07  2.07E-06
 
 OM24
+        2.55E-06 -4.02E-07 -6.55E-07  7.11E-08  5.02E-06 -2.18E-06 -1.25E-06  9.55E-07 -4.55E-08  8.11E-07 -3.64E-07  1.97E-07
          3.57E-07  4.64E-07  5.26E-07  5.21E-07  1.13E-06  1.34E-06
 
 OM33
+        4.38E-06  1.66E-06 -1.33E-06 -4.45E-07  1.01E-05 -5.26E-06 -2.94E-06  1.18E-06 -1.27E-06 -1.11E-06 -3.17E-06  4.44E-07
          4.60E-07  1.29E-06  6.62E-07  2.92E-07  1.53E-06  7.05E-07  7.99E-06
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        3.02E-06  2.74E-06 -9.63E-07 -7.27E-07  7.77E-06  2.18E-06 -2.43E-06 -5.36E-07 -2.93E-07  4.28E-07 -2.29E-06  4.18E-07
          3.94E-07  9.03E-07  8.47E-07  2.28E-07  1.21E-06  7.91E-07  4.54E-06  3.94E-06
 
 OM44
+       -1.57E-07  3.53E-06 -2.41E-08 -7.31E-07  8.55E-06  5.78E-06 -2.63E-06 -1.04E-06  1.18E-06  3.19E-06 -1.83E-06  4.13E-07
          3.77E-07  7.35E-07  9.64E-07  9.82E-08  8.20E-07  9.46E-07  2.41E-06  3.05E-06  4.00E-06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (S)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        3.80E-02
 
 TH 2
+        2.84E-03  2.88E-02
 
 TH 3
+       -9.77E-01  4.01E-03  1.11E-02
 
 TH 4
+       -1.34E-03 -9.55E-01 -1.32E-03  8.42E-03
 
 TH 5
+        3.49E-01  3.68E-02 -3.11E-01 -1.03E-02  4.84E-02
 
 TH 6
+       -1.20E-02  1.80E-01  1.37E-02 -1.51E-01  1.13E-02  4.10E-02
 
 TH 7
+       -3.28E-01 -3.96E-02  2.95E-01  1.41E-02 -9.79E-01 -1.78E-02  1.34E-02
 
 TH 8
+        1.36E-02 -1.54E-01 -1.24E-02  1.48E-01  7.16E-03 -9.68E-01 -6.28E-04  1.19E-02
 
 TH 9
+        4.37E-02  1.98E-01 -7.15E-04 -1.03E-01  2.27E-01  1.17E-01 -2.30E-01 -5.65E-02  1.05E-02
 
 TH10
+        4.06E-02  1.43E-01 -6.31E-04 -6.35E-02  2.49E-01  1.31E-01 -2.47E-01 -5.12E-02  7.01E-01  8.94E-03
 
 TH11
+       -6.68E-02 -4.24E-02  8.10E-02  5.12E-02 -6.05E-02  6.06E-03  6.22E-02 -2.26E-03  5.50E-02  5.48E-02  2.98E-03
 
 OM11
+        1.08E-01 -5.32E-03 -1.19E-01 -2.04E-02  2.62E-02  6.06E-03 -3.39E-02 -3.49E-03 -3.75E-02  1.37E-02 -1.42E-01  1.01E-03
 
 OM12
+        4.68E-02  6.53E-03 -6.95E-02 -7.12E-03  4.52E-02  2.09E-02 -3.23E-02  1.37E-03 -1.77E-02  4.34E-02 -2.18E-01  4.23E-01
          8.65E-04
 
 OM13
+        9.53E-03  4.92E-02 -2.94E-02 -5.34E-02  3.20E-02  4.55E-04 -3.24E-02  3.60E-03  5.45E-02  5.69E-02 -1.34E-01  4.10E-01
          2.97E-01  1.29E-03
 
 OM14
+        4.56E-02 -1.43E-02 -6.00E-02  3.14E-02  1.12E-01 -4.11E-03 -1.01E-01  1.04E-02  6.43E-02  9.23E-02 -1.02E-01  3.84E-01
          3.06E-01  7.15E-01  1.12E-03
 
 OM22
+        6.70E-02  1.56E-02 -6.45E-02 -7.22E-03  7.45E-02 -1.16E-01 -3.64E-02  1.20E-01 -2.28E-03 -3.43E-03 -2.41E-01  6.04E-02
          3.60E-01  3.32E-02  8.92E-02  1.50E-03
 
 OM23
+        5.66E-03  1.35E-02 -1.47E-02 -2.36E-02 -6.50E-03  2.13E-02  9.61E-03 -3.82E-03 -4.62E-02 -2.76E-02 -1.60E-01  2.06E-01
          3.43E-01  4.22E-01  3.14E-01  2.46E-01  1.44E-03
 
 OM24
+        5.81E-02 -1.21E-02 -5.11E-02  7.30E-03  8.97E-02 -4.59E-02 -8.08E-02  6.97E-02 -3.74E-03  7.85E-02 -1.06E-01  1.69E-01
          3.57E-01  3.11E-01  4.07E-01  3.01E-01  6.77E-01  1.16E-03
 
 OM33
+        4.07E-02  2.03E-02 -4.24E-02 -1.87E-02  7.38E-02 -4.53E-02 -7.75E-02  3.51E-02 -4.28E-02 -4.38E-02 -3.76E-01  1.56E-01
          1.88E-01  3.52E-01  2.10E-01  6.90E-02  3.76E-01  2.16E-01  2.83E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        4.00E-02  4.79E-02 -4.37E-02 -4.35E-02  8.09E-02  2.67E-02 -9.11E-02 -2.28E-02 -1.40E-02  2.41E-02 -3.87E-01  2.09E-01
          2.30E-01  3.52E-01  3.82E-01  7.66E-02  4.22E-01  3.45E-01  8.09E-01  1.98E-03
 
 OM44
+       -2.06E-03  6.12E-02 -1.09E-03 -4.34E-02  8.84E-02  7.05E-02 -9.79E-02 -4.40E-02  5.59E-02  1.79E-01 -3.07E-01  2.05E-01
          2.18E-01  2.85E-01  4.31E-01  3.28E-02  2.85E-01  4.09E-01  4.26E-01  7.68E-01  2.00E-03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          ITERATIVE TWO STAGE (NO PRIOR)                        ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (S)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.72E+04
 
 TH 2
+        0.00E+00  1.71E+04
 
 TH 3
+        5.67E+04  0.00E+00  1.97E+05
 
 TH 4
+        0.00E+00  5.55E+04  0.00E+00  1.95E+05
 
 TH 5
+       -1.95E+03  0.00E+00 -5.23E+03  0.00E+00  1.15E+04
 
 TH 6
+        0.00E+00 -3.68E+03  0.00E+00 -1.19E+04  0.00E+00  1.16E+04
 
 TH 7
+       -5.23E+03  0.00E+00 -1.54E+04  0.00E+00  4.01E+04  0.00E+00  1.46E+05
 
 TH 8
+        0.00E+00 -1.19E+04  0.00E+00 -3.97E+04  0.00E+00  3.83E+04  0.00E+00  1.34E+05
 
 TH 9
+       -1.68E+03 -3.54E+03 -5.81E+03 -1.01E+04  5.58E+02  8.82E+01  2.81E+03  4.46E+02  1.95E+04
 
 TH10
+       -1.08E+03 -1.20E+02 -4.07E+03 -7.26E+02 -6.85E+02 -3.29E+03 -8.78E+02 -1.07E+04 -1.49E+04  2.82E+04
 
 TH11
+       -1.98E+03 -1.31E+03 -8.55E+03 -5.86E+03 -1.08E+03  2.02E+02 -4.93E+03  1.34E+02  2.77E+02 -3.90E+03  1.49E+05
 
 OM11
+       -1.02E+04  1.38E+04 -2.14E+04  4.69E+04  9.97E+03 -4.18E+02  3.30E+04 -2.50E+03  7.71E+03 -4.81E+03  1.69E+04  1.43E+06
 
 OM12
+        2.09E+04 -5.33E+02  6.93E+04 -4.18E+03 -5.30E+03 -1.17E+04 -1.86E+04 -3.27E+04  5.64E+03 -9.05E+03  3.84E+04 -5.88E+05
          2.10E+06
 
 OM13
+        7.63E+03  9.73E+03  2.38E+04  4.88E+04  6.19E+03 -1.44E+03  1.68E+04 -6.39E+03 -8.78E+03 -4.46E+03  5.99E+03 -2.42E+05
         -3.84E+04  1.67E+06
 
 OM14
+        8.90E+03 -1.86E+04  2.94E+04 -8.25E+04 -1.83E+04  1.02E+04 -5.46E+04  3.36E+04 -4.97E+03  8.32E+03 -3.12E+04 -2.01E+05
         -6.71E+04 -1.30E+06  2.29E+06
 
 OM22
+       -1.85E+03 -5.83E+03 -5.91E+03 -1.63E+04 -1.38E+04  4.61E+03 -4.82E+04  7.37E+03 -2.49E+03  8.51E+02  6.94E+04  5.92E+04
         -3.19E+05  9.07E+04 -2.92E+04  6.17E+05
 
 OM23
+        6.93E+03 -2.29E+01  1.94E+04 -1.44E+03 -1.05E+02 -8.77E+03 -5.56E+03 -2.13E+04 -1.21E+03  9.31E+03 -2.18E+03 -2.51E+04
         -8.61E+04 -4.06E+05  3.12E+05 -4.71E+04  1.21E+06
 
 OM24
+       -1.79E+04  1.13E+04 -5.37E+04  3.72E+04  3.36E+03 -6.48E+02  1.70E+04 -1.34E+04  8.96E+03 -1.06E+04 -4.72E+04  1.34E+05
         -2.14E+05  2.46E+05 -4.86E+05 -1.73E+05 -9.39E+05  1.89E+06
 
 OM33
+       -1.13E+03 -6.99E+03 -6.38E+03 -2.75E+04 -5.53E+03  5.31E+03 -1.61E+04  1.22E+04  3.54E+03 -4.44E+03  4.14E+04  2.03E+04
         -4.01E+04 -2.79E+05  2.98E+05  3.47E+04  3.09E+04 -4.76E+04  5.71E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+       -6.38E+03  8.79E+03 -1.31E+04  3.62E+04  1.13E+04 -9.60E+02  3.87E+04  3.11E+03 -9.02E+03  2.46E+04  4.84E+02 -1.97E+03
          5.10E+04  3.41E+05 -4.61E+05 -5.23E+04 -3.71E+05  3.15E+05 -9.05E+05  2.24E+06
 
 OM44
+        6.10E+03 -5.88E+03  1.24E+04 -1.84E+04 -1.54E+03 -4.62E+03 -5.36E+03 -1.46E+04  1.25E+04 -3.25E+04  5.59E+04 -3.25E+04
         -2.07E+04 -2.45E+04 -8.20E+04  9.57E+04  2.45E+05 -4.01E+05  3.50E+05 -1.13E+06  1.03E+06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 
 
 #TBLN:      2
 #METH: Importance Sampling (No Prior)
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            2208
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      8
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     8
 NOPRIOR SETTING (NOPRIOR):                 ON
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example2.ext
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
 GRADIENT/GIBBS PATTERN (GRD):              DDDDDDDDDDS
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

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8   9  10
 THETAS THAT ARE SIGMA-LIKE:
  11
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=  -10780.2980251398 eff.=     304. Smpl.=     300. Fit.= 0.98310
 iteration            1 OBJ=  -10780.6610290265 eff.=     124. Smpl.=     300. Fit.= 0.90173
 iteration            2 OBJ=  -10780.8233352404 eff.=     119. Smpl.=     300. Fit.= 0.89854
 iteration            3 OBJ=  -10782.4559669737 eff.=     119. Smpl.=     300. Fit.= 0.89841
 iteration            4 OBJ=  -10780.4700626909 eff.=     122. Smpl.=     300. Fit.= 0.90070
 iteration            5 OBJ=  -10782.2108096674 eff.=     120. Smpl.=     300. Fit.= 0.89942
 iteration            6 OBJ=  -10783.5226777336 eff.=     121. Smpl.=     300. Fit.= 0.89937
 iteration            7 OBJ=  -10781.2203006768 eff.=     121. Smpl.=     300. Fit.= 0.89996
 iteration            8 OBJ=  -10781.9272924653 eff.=     121. Smpl.=     300. Fit.= 0.89996
 iteration            9 OBJ=  -10782.5370946805 eff.=     120. Smpl.=     300. Fit.= 0.89886
 iteration           10 OBJ=  -10784.3160132351 eff.=     121. Smpl.=     300. Fit.= 0.89964
 iteration           11 OBJ=  -10781.8799403199 eff.=     122. Smpl.=     300. Fit.= 0.90113
 iteration           12 OBJ=  -10783.6229290287 eff.=     120. Smpl.=     300. Fit.= 0.89909
 iteration           13 OBJ=  -10783.3543023934 eff.=     120. Smpl.=     300. Fit.= 0.89949
 iteration           14 OBJ=  -10782.8263529080 eff.=     120. Smpl.=     300. Fit.= 0.89877
 iteration           15 OBJ=  -10782.1900732069 eff.=     121. Smpl.=     300. Fit.= 0.90002
 iteration           16 OBJ=  -10780.4115294156 eff.=     121. Smpl.=     300. Fit.= 0.90057
 iteration           17 OBJ=  -10781.8699691324 eff.=     119. Smpl.=     300. Fit.= 0.89821
 iteration           18 OBJ=  -10780.3233988869 eff.=     122. Smpl.=     300. Fit.= 0.90080
 iteration           19 OBJ=  -10779.9849555581 eff.=     121. Smpl.=     300. Fit.= 0.90022
 iteration           20 OBJ=  -10782.1904878407 eff.=     120. Smpl.=     300. Fit.= 0.89933
 iteration           21 OBJ=  -10783.1058987664 eff.=     120. Smpl.=     300. Fit.= 0.89878
 iteration           22 OBJ=  -10783.3814576766 eff.=     124. Smpl.=     300. Fit.= 0.90174
 iteration           23 OBJ=  -10779.1113000283 eff.=     118. Smpl.=     300. Fit.= 0.89853
 iteration           24 OBJ=  -10781.0946734138 eff.=     122. Smpl.=     300. Fit.= 0.90095
 iteration           25 OBJ=  -10779.8247696467 eff.=     120. Smpl.=     300. Fit.= 0.89951
 iteration           26 OBJ=  -10779.9803211255 eff.=     120. Smpl.=     300. Fit.= 0.90000
 iteration           27 OBJ=  -10785.9629809746 eff.=     122. Smpl.=     300. Fit.= 0.89974
 Convergence achieved
 iteration           27 OBJ=  -10779.7809930198 eff.=     121. Smpl.=     300. Fit.= 0.90054
 
 #TERM:
 OPTIMIZATION WAS COMPLETED


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -1.0108E-04 -1.0426E-04  5.6875E-05 -1.1009E-04
 SE:             4.7137E-03  2.9745E-03  2.8988E-03  3.6761E-03
 N:                     400         400         400         400
 
 P VAL.:         9.8289E-01  9.7204E-01  9.8435E-01  9.7611E-01
 
 ETASHRINKSD(%)  7.0444E+00  3.3496E+01  4.1481E+01  2.5326E+01
 ETASHRINKVR(%)  1.3593E+01  5.5773E+01  6.5755E+01  4.4238E+01
 EBVSHRINKSD(%)  7.0448E+00  3.3554E+01  4.1505E+01  2.5196E+01
 EBVSHRINKVR(%)  1.3593E+01  5.5850E+01  6.5783E+01  4.4044E+01
 EPSSHRINKSD(%)  2.6387E+01
 EPSSHRINKVR(%)  4.5811E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2000
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    3675.75413281869     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -10779.7809930198     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -7104.02686020115     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1600
  
 #TERE:
 Elapsed estimation  time in seconds:    55.24
 Elapsed covariance  time in seconds:     5.32
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          IMPORTANCE SAMPLING (NO PRIOR)                        ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -10779.781       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          IMPORTANCE SAMPLING (NO PRIOR)                        ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.30E+00  3.25E+00 -6.12E-01 -2.08E-01  7.32E-01  1.14E+00  3.36E-01  1.92E-01  6.90E-01  2.30E+00  1.00E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.03E-02
 
 ETA2
+        1.85E-04  8.02E-03
 
 ETA3
+        1.21E-03 -1.82E-04  9.84E-03
 
 ETA4
+       -6.41E-04  5.03E-04  1.88E-03  9.72E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.02E-01
 
 ETA2
+        2.03E-02  8.96E-02
 
 ETA3
+        1.20E-01 -2.05E-02  9.92E-02
 
 ETA4
+       -6.41E-02  5.70E-02  1.92E-01  9.86E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          IMPORTANCE SAMPLING (NO PRIOR)                        ********************
 ********************                          STANDARD ERROR OF ESTIMATE (R)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.27E-02  2.88E-02  9.54E-03  8.36E-03  3.95E-02  3.65E-02  1.14E-02  1.05E-02  1.05E-02  8.62E-03  2.80E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        9.73E-04
 
 ETA2
+        8.30E-04  1.36E-03
 
 ETA3
+        1.31E-03  1.40E-03  3.31E-03
 
 ETA4
+        1.02E-03  1.08E-03  2.29E-03  1.96E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        4.79E-03
 
 ETA2
+        9.05E-02  7.57E-03
 
 ETA3
+        1.19E-01  1.59E-01  1.67E-02
 
 ETA4
+        1.06E-01  1.20E-01  1.91E-01  9.95E-03
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          IMPORTANCE SAMPLING (NO PRIOR)                        ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (R)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.07E-03
 
 TH 2
+        1.56E-05  8.31E-04
 
 TH 3
+       -3.02E-04 -2.89E-06  9.11E-05
 
 TH 4
+       -1.03E-06 -2.31E-04  1.09E-07  6.99E-05
 
 TH 5
+        2.94E-04  2.62E-05 -7.92E-05 -2.71E-06  1.56E-03
 
 TH 6
+        3.01E-05  3.11E-04 -6.02E-06 -8.12E-05  4.49E-05  1.33E-03
 
 TH 7
+       -8.18E-05 -8.37E-06  2.30E-05  9.12E-07 -4.36E-04 -1.41E-05  1.30E-04
 
 TH 8
+       -5.03E-06 -8.33E-05  1.04E-06  2.38E-05 -6.59E-06 -3.69E-04  2.18E-06  1.11E-04
 
 TH 9
+        4.35E-05  4.60E-05 -9.21E-06 -4.78E-06  5.89E-05  7.36E-05 -1.86E-05 -1.24E-05  1.10E-04
 
 TH10
+        2.71E-05  3.21E-05 -5.02E-06 -3.64E-06  5.32E-05  6.16E-05 -1.56E-05 -9.58E-06  6.34E-05  7.43E-05
 
 TH11
+        1.06E-06 -3.19E-07 -3.31E-07  1.29E-07  3.34E-06 -2.73E-06 -1.02E-06  8.35E-07  1.27E-06  9.18E-07  7.83E-06
 
 OM11
+        3.62E-07 -1.98E-07 -3.13E-08  8.92E-08  4.75E-07  4.13E-07 -3.69E-08 -7.72E-09  3.80E-07  3.33E-07 -4.14E-07  9.47E-07
 
 OM12
+        4.44E-08 -3.12E-07  5.17E-08  1.49E-07 -2.78E-08 -3.68E-07  1.58E-07  2.29E-07  3.01E-07  2.75E-07 -3.93E-07  2.79E-07
          6.89E-07
 
 OM13
+        9.29E-07 -9.72E-08 -7.97E-08  8.83E-08  2.11E-06  1.23E-06 -4.42E-07 -1.06E-07  8.20E-07  8.22E-07 -8.28E-07  5.60E-07
          2.25E-07  1.73E-06
 
 OM14
+        1.21E-06 -1.56E-06 -2.72E-07  5.47E-07  9.81E-07  8.12E-07 -1.76E-07 -8.14E-08  6.38E-07  4.31E-07 -6.38E-07  3.71E-07
          2.13E-07  9.49E-07  1.04E-06
 
 OM22
+       -7.37E-07 -4.12E-07  2.58E-07  1.19E-07 -2.23E-06 -1.75E-06  8.65E-07  5.62E-07 -5.95E-07 -4.12E-07 -1.08E-06  7.72E-08
          3.47E-07  5.60E-08  5.32E-08  1.84E-06
 
 OM23
+        8.79E-07 -1.49E-08 -1.64E-07  8.61E-08  2.62E-06 -9.93E-07 -4.69E-07  4.20E-07  1.12E-06  7.63E-07 -5.03E-07  1.81E-07
          3.31E-07  6.71E-07  3.42E-07  2.99E-07  1.97E-06
 
 OM24
+        2.09E-07 -4.85E-07 -1.65E-08  1.61E-07  1.39E-06 -2.97E-06 -2.66E-07  9.33E-07  4.33E-07  2.94E-07 -3.04E-07  9.62E-08
          2.15E-07  2.61E-07  3.00E-07  3.72E-07  9.18E-07  1.18E-06
 
 OM33
+        1.57E-06  1.69E-06 -1.81E-07 -1.66E-07  4.14E-06  6.42E-06 -9.75E-07 -1.25E-06  1.08E-06  1.22E-06 -3.85E-06  5.34E-07
          2.98E-07  2.31E-06  1.24E-06  2.11E-07  1.61E-06  4.74E-07  1.10E-05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        8.31E-07  1.87E-06 -7.80E-08 -3.44E-07  1.90E-06  6.07E-06 -4.21E-07 -1.40E-06  7.35E-07  5.55E-07 -2.67E-06  3.45E-07
          2.04E-07  1.44E-06  1.01E-06  8.82E-08  1.10E-06  4.67E-07  6.63E-06  5.23E-06
 
 OM44
+        1.46E-07  1.69E-06  4.35E-08 -3.90E-07  1.22E-07  5.52E-06  3.12E-08 -1.42E-06  1.47E-07  1.96E-07 -1.95E-06  2.15E-07
          1.32E-07  8.15E-07  7.69E-07  4.64E-08  5.73E-07  5.15E-07  3.70E-06  3.62E-06  3.85E-06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          IMPORTANCE SAMPLING (NO PRIOR)                        ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (R)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        3.27E-02
 
 TH 2
+        1.65E-02  2.88E-02
 
 TH 3
+       -9.68E-01 -1.05E-02  9.54E-03
 
 TH 4
+       -3.75E-03 -9.59E-01  1.37E-03  8.36E-03
 
 TH 5
+        2.28E-01  2.30E-02 -2.10E-01 -8.21E-03  3.95E-02
 
 TH 6
+        2.52E-02  2.96E-01 -1.73E-02 -2.66E-01  3.12E-02  3.65E-02
 
 TH 7
+       -2.19E-01 -2.55E-02  2.11E-01  9.58E-03 -9.69E-01 -3.38E-02  1.14E-02
 
 TH 8
+       -1.46E-02 -2.74E-01  1.03E-02  2.70E-01 -1.58E-02 -9.59E-01  1.81E-02  1.05E-02
 
 TH 9
+        1.26E-01  1.52E-01 -9.18E-02 -5.45E-02  1.42E-01  1.92E-01 -1.56E-01 -1.12E-01  1.05E-02
 
 TH10
+        9.60E-02  1.29E-01 -6.11E-02 -5.05E-02  1.56E-01  1.96E-01 -1.58E-01 -1.05E-01  7.00E-01  8.62E-03
 
 TH11
+        1.15E-02 -3.95E-03 -1.24E-02  5.52E-03  3.02E-02 -2.68E-02 -3.19E-02  2.83E-02  4.33E-02  3.81E-02  2.80E-03
 
 OM11
+        1.14E-02 -7.05E-03 -3.37E-03  1.10E-02  1.24E-02  1.16E-02 -3.33E-03 -7.52E-04  3.72E-02  3.97E-02 -1.52E-01  9.73E-04
 
 OM12
+        1.63E-03 -1.30E-02  6.52E-03  2.15E-02 -8.50E-04 -1.22E-02  1.67E-02  2.61E-02  3.45E-02  3.84E-02 -1.69E-01  3.45E-01
          8.30E-04
 
 OM13
+        2.16E-02 -2.57E-03 -6.35E-03  8.04E-03  4.06E-02  2.57E-02 -2.95E-02 -7.64E-03  5.94E-02  7.26E-02 -2.25E-01  4.38E-01
          2.06E-01  1.31E-03
 
 OM14
+        3.62E-02 -5.30E-02 -2.80E-02  6.42E-02  2.44E-02  2.18E-02 -1.52E-02 -7.58E-03  5.96E-02  4.91E-02 -2.24E-01  3.74E-01
          2.52E-01  7.09E-01  1.02E-03
 
 OM22
+       -1.66E-02 -1.05E-02  1.99E-02  1.05E-02 -4.16E-02 -3.54E-02  5.60E-02  3.92E-02 -4.18E-02 -3.53E-02 -2.85E-01  5.85E-02
          3.08E-01  3.14E-02  3.85E-02  1.36E-03
 
 OM23
+        1.91E-02 -3.69E-04 -1.23E-02  7.34E-03  4.73E-02 -1.94E-02 -2.93E-02  2.83E-02  7.63E-02  6.30E-02 -1.28E-01  1.32E-01
          2.84E-01  3.64E-01  2.39E-01  1.57E-01  1.40E-03
 
 OM24
+        5.89E-03 -1.55E-02 -1.59E-03  1.77E-02  3.25E-02 -7.49E-02 -2.15E-02  8.16E-02  3.80E-02  3.14E-02 -1.00E-01  9.11E-02
          2.38E-01  1.83E-01  2.72E-01  2.53E-01  6.03E-01  1.08E-03
 
 OM33
+        1.45E-02  1.77E-02 -5.73E-03 -6.01E-03  3.17E-02  5.31E-02 -2.58E-02 -3.57E-02  3.11E-02  4.26E-02 -4.15E-01  1.66E-01
          1.08E-01  5.30E-01  3.68E-01  4.68E-02  3.46E-01  1.32E-01  3.31E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        1.11E-02  2.84E-02 -3.57E-03 -1.80E-02  2.11E-02  7.28E-02 -1.62E-02 -5.82E-02  3.06E-02  2.81E-02 -4.18E-01  1.55E-01
          1.07E-01  4.78E-01  4.33E-01  2.84E-02  3.43E-01  1.88E-01  8.75E-01  2.29E-03
 
 OM44
+        2.28E-03  2.99E-02  2.32E-03 -2.38E-02  1.57E-03  7.70E-02  1.39E-03 -6.87E-02  7.12E-03  1.16E-02 -3.54E-01  1.13E-01
          8.13E-02  3.16E-01  3.85E-01  1.74E-02  2.08E-01  2.42E-01  5.68E-01  8.07E-01  1.96E-03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                          IMPORTANCE SAMPLING (NO PRIOR)                        ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (R)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.58E+04
 
 TH 2
+        0.00E+00  1.84E+04
 
 TH 3
+        5.20E+04  0.00E+00  1.83E+05
 
 TH 4
+        0.00E+00  6.03E+04  0.00E+00  2.13E+05
 
 TH 5
+       -2.11E+03  0.00E+00 -6.69E+03  0.00E+00  1.11E+04
 
 TH 6
+        0.00E+00 -3.51E+03  0.00E+00 -1.13E+04  0.00E+00  1.15E+04
 
 TH 7
+       -6.69E+03  0.00E+00 -2.31E+04  0.00E+00  3.73E+04  0.00E+00  1.33E+05
 
 TH 8
+        0.00E+00 -1.13E+04  0.00E+00 -3.95E+04  0.00E+00  3.75E+04  0.00E+00  1.33E+05
 
 TH 9
+       -1.31E+03 -3.83E+03 -3.49E+03 -1.23E+04  1.10E+03 -6.20E+02  4.25E+03 -1.41E+03  1.92E+04
 
 TH10
+       -9.64E+02 -2.94E+02 -3.52E+03 -9.05E+02 -6.96E+02 -3.15E+03 -1.31E+03 -9.71E+03 -1.46E+04  2.78E+04
 
 TH11
+       -2.65E+02 -2.68E+02 -7.52E+02 -1.18E+03 -1.05E+03 -8.89E+02 -2.94E+03 -3.67E+03 -1.05E+03 -1.08E+03  1.77E+05
 
 OM11
+        5.23E+01  4.78E+02  4.21E+02  1.72E+03 -1.26E+03 -5.65E+02 -4.56E+03 -1.58E+03 -1.54E+03  3.76E+02  3.41E+04  1.46E+06
 
 OM12
+       -1.92E+03 -1.39E+03 -7.92E+03 -4.71E+03 -3.48E+03 -3.30E+03 -1.32E+04 -1.26E+04  2.30E+02 -3.46E+03  2.26E+04 -4.82E+05
          1.96E+06
 
 OM13
+       -6.05E+03  5.41E+03 -2.39E+04  2.68E+04 -6.73E+01 -2.27E+03  1.43E+03 -1.04E+04  3.93E+03 -7.05E+03 -2.39E+04 -4.37E+05
          1.58E+05  1.70E+06
 
 OM14
+        3.28E+03 -3.34E+03  1.74E+04 -2.99E+04 -2.01E+03  3.21E+02 -8.50E+03  5.29E+03 -7.39E+03  5.74E+03  3.17E+04 -8.70E+04
         -3.00E+05 -1.26E+06  2.38E+06
 
 OM22
+       -9.06E+02 -1.14E+03 -3.16E+03 -3.93E+03 -3.01E+03 -1.68E+03 -1.34E+04 -6.26E+03  2.31E+03  1.82E+03  9.80E+04  4.82E+04
         -2.92E+05 -1.31E+04  8.19E+04  6.92E+05
 
 OM23
+        2.22E+03 -1.29E+02  8.93E+03 -2.74E+03 -5.42E+03  2.31E+03 -1.84E+04  9.83E+03 -5.97E+03  7.47E+02 -1.54E+04  7.67E+04
         -2.63E+05 -3.46E+05  3.20E+05  3.11E+04  1.05E+06
 
 OM24
+       -7.60E+02  1.91E+03 -4.62E+03  1.05E+04  1.39E+03  8.33E+02  7.86E+03 -8.67E+03 -4.68E+01 -1.92E+03 -1.91E+04  6.72E+02
          1.38E+04  2.56E+05 -4.62E+05 -2.12E+05 -7.78E+05  1.62E+06
 
 OM33
+       -1.89E-01 -1.75E+03  6.51E+02 -8.27E+03 -3.91E+02 -1.04E+03 -5.45E+02 -4.07E+03  3.37E+03 -4.41E+03  4.65E+04  3.37E+04
         -9.37E+03 -2.82E+05  2.21E+05 -2.70E+02  1.32E+04 -1.14E+04  5.68E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+       -5.60E+02  5.49E+02 -2.92E+03  5.94E+03  9.91E+02 -9.74E+02  3.68E+03 -2.64E+03 -6.79E+03  9.50E+03  5.35E+03  1.36E+04
          5.86E+04  1.88E+05 -3.60E+05  1.77E+04 -2.77E+05  2.48E+05 -8.97E+05  2.17E+06
 
 OM44
+       -1.05E+02 -5.92E+02 -4.33E+02 -1.61E+03 -3.29E+02 -4.11E+02 -2.31E+03  2.70E+03  5.04E+03 -4.46E+03  3.98E+04  4.94E+03
         -7.43E+03  1.09E+04 -3.80E+04  4.26E+04  2.02E+05 -2.97E+05  3.35E+05 -1.14E+06  1.05E+06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 
 
 #TBLN:      3
 #METH: Stochastic Approximation Expectation-Maximization (No Prior)
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            2208
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      8
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     8
 NOPRIOR SETTING (NOPRIOR):                 ON
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example2.ext
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
 GRADIENT/GIBBS PATTERN (GRD):              DDDDDDDDDDS
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  3
 CONVERGENCE INTERVAL (CINTERVAL):          10
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 BURN-IN ITERATIONS (NBURN):                3000
 ITERATIONS (NITER):                        2000
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       11456
 MC SAMPLES PER SUBJECT (ISAMPLE):          2
 RANDOM SAMPLING METHOD (RANMETHOD):        3U
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

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8   9  10
 THETAS THAT ARE SIGMA-LIKE:
  11
 
 MONITORING OF SEARCH:

 Stochastic/Burn-in Mode
 iteration        -3000 SAEMOBJ=  -19231.0174739001
 iteration        -2990 SAEMOBJ=  -19387.3278561200
 iteration        -2980 SAEMOBJ=  -19418.5417122258
 iteration        -2970 SAEMOBJ=  -19445.7492871097
 iteration        -2960 SAEMOBJ=  -19437.6170595463
 iteration        -2950 SAEMOBJ=  -19497.0747751477
 iteration        -2940 SAEMOBJ=  -19467.1628605598
 iteration        -2930 SAEMOBJ=  -19477.1704586428
 iteration        -2920 SAEMOBJ=  -19463.7467651272
 iteration        -2910 SAEMOBJ=  -19468.6557892722
 iteration        -2900 SAEMOBJ=  -19439.9861923145
 iteration        -2890 SAEMOBJ=  -19619.9547859568
 Convergence achieved
 Reduced Stochastic/Accumulation Mode
 iteration            0 SAEMOBJ=  -19557.5325885114
 iteration           10 SAEMOBJ=  -19700.0360719512
 iteration           20 SAEMOBJ=  -19701.6840378869
 iteration           30 SAEMOBJ=  -19697.1399168863
 iteration           40 SAEMOBJ=  -19694.1854625138
 iteration           50 SAEMOBJ=  -19691.4850408402
 iteration           60 SAEMOBJ=  -19688.3081064560
 iteration           70 SAEMOBJ=  -19687.1970594812
 iteration           80 SAEMOBJ=  -19686.4903082082
 iteration           90 SAEMOBJ=  -19686.4293523707
 iteration          100 SAEMOBJ=  -19685.8267415077
 iteration          110 SAEMOBJ=  -19684.6510269996
 iteration          120 SAEMOBJ=  -19683.6358703418
 iteration          130 SAEMOBJ=  -19683.8596538376
 iteration          140 SAEMOBJ=  -19684.3466087238
 iteration          150 SAEMOBJ=  -19684.8402670651
 iteration          160 SAEMOBJ=  -19685.3031639363
 iteration          170 SAEMOBJ=  -19685.1851599943
 iteration          180 SAEMOBJ=  -19685.4353348545
 iteration          190 SAEMOBJ=  -19685.3989172726
 iteration          200 SAEMOBJ=  -19685.6074354025
 iteration          210 SAEMOBJ=  -19685.9376290948
 iteration          220 SAEMOBJ=  -19685.8577144725
 iteration          230 SAEMOBJ=  -19686.2045607318
 iteration          240 SAEMOBJ=  -19686.6812736439
 iteration          250 SAEMOBJ=  -19687.1901077744
 iteration          260 SAEMOBJ=  -19686.8996251235
 iteration          270 SAEMOBJ=  -19687.2215324757
 iteration          280 SAEMOBJ=  -19687.6992131768
 iteration          290 SAEMOBJ=  -19688.0724082399
 iteration          300 SAEMOBJ=  -19688.4905317777
 iteration          310 SAEMOBJ=  -19688.7531590835
 iteration          320 SAEMOBJ=  -19688.7210939407
 iteration          330 SAEMOBJ=  -19689.0666501710
 iteration          340 SAEMOBJ=  -19689.2069823266
 iteration          350 SAEMOBJ=  -19689.2144645008
 iteration          360 SAEMOBJ=  -19689.1390145867
 iteration          370 SAEMOBJ=  -19689.0801581917
 iteration          380 SAEMOBJ=  -19689.2319126503
 iteration          390 SAEMOBJ=  -19689.0071968747
 iteration          400 SAEMOBJ=  -19688.7777503334
 iteration          410 SAEMOBJ=  -19689.0606530074
 iteration          420 SAEMOBJ=  -19689.0354049805
 iteration          430 SAEMOBJ=  -19688.9821973313
 iteration          440 SAEMOBJ=  -19688.7029248225
 iteration          450 SAEMOBJ=  -19688.4934547207
 iteration          460 SAEMOBJ=  -19688.4262902909
 iteration          470 SAEMOBJ=  -19688.4480046876
 iteration          480 SAEMOBJ=  -19688.3347490744
 iteration          490 SAEMOBJ=  -19688.2799157102
 iteration          500 SAEMOBJ=  -19688.6063440868
 iteration          510 SAEMOBJ=  -19688.4377200515
 iteration          520 SAEMOBJ=  -19688.3903977631
 iteration          530 SAEMOBJ=  -19688.1371081096
 iteration          540 SAEMOBJ=  -19687.8846933461
 iteration          550 SAEMOBJ=  -19687.9072496530
 iteration          560 SAEMOBJ=  -19687.8622408996
 iteration          570 SAEMOBJ=  -19687.8250995603
 iteration          580 SAEMOBJ=  -19687.8941894944
 iteration          590 SAEMOBJ=  -19687.9164143527
 iteration          600 SAEMOBJ=  -19687.9589959852
 iteration          610 SAEMOBJ=  -19688.1685264472
 iteration          620 SAEMOBJ=  -19688.1294839307
 iteration          630 SAEMOBJ=  -19688.0905593074
 iteration          640 SAEMOBJ=  -19688.2211571643
 iteration          650 SAEMOBJ=  -19688.1771292468
 iteration          660 SAEMOBJ=  -19688.2940487110
 iteration          670 SAEMOBJ=  -19688.2861641664
 iteration          680 SAEMOBJ=  -19688.3030501822
 iteration          690 SAEMOBJ=  -19688.4030676318
 iteration          700 SAEMOBJ=  -19688.4264851816
 iteration          710 SAEMOBJ=  -19688.4653915549
 iteration          720 SAEMOBJ=  -19688.5511262271
 iteration          730 SAEMOBJ=  -19688.4619339304
 iteration          740 SAEMOBJ=  -19688.3662147660
 iteration          750 SAEMOBJ=  -19688.5404059089
 iteration          760 SAEMOBJ=  -19688.3461543379
 iteration          770 SAEMOBJ=  -19688.1765680731
 iteration          780 SAEMOBJ=  -19688.2017134958
 iteration          790 SAEMOBJ=  -19688.1113919829
 iteration          800 SAEMOBJ=  -19688.1660169192
 iteration          810 SAEMOBJ=  -19688.3766571465
 iteration          820 SAEMOBJ=  -19688.4774174446
 iteration          830 SAEMOBJ=  -19688.4737579437
 iteration          840 SAEMOBJ=  -19688.6179048974
 iteration          850 SAEMOBJ=  -19688.4586399634
 iteration          860 SAEMOBJ=  -19688.4074998923
 iteration          870 SAEMOBJ=  -19688.4272772840
 iteration          880 SAEMOBJ=  -19688.2622910593
 iteration          890 SAEMOBJ=  -19688.4365352365
 iteration          900 SAEMOBJ=  -19688.4999129910
 iteration          910 SAEMOBJ=  -19688.5214901865
 iteration          920 SAEMOBJ=  -19688.4582839471
 iteration          930 SAEMOBJ=  -19688.4200327429
 iteration          940 SAEMOBJ=  -19688.3527181914
 iteration          950 SAEMOBJ=  -19688.3274139807
 iteration          960 SAEMOBJ=  -19688.2554907329
 iteration          970 SAEMOBJ=  -19688.2827295293
 iteration          980 SAEMOBJ=  -19688.2491377553
 iteration          990 SAEMOBJ=  -19688.3395946278
 iteration         1000 SAEMOBJ=  -19688.2768007275
 iteration         1010 SAEMOBJ=  -19688.2094223278
 iteration         1020 SAEMOBJ=  -19688.1459084371
 iteration         1030 SAEMOBJ=  -19688.3597995622
 iteration         1040 SAEMOBJ=  -19688.2866981688
 iteration         1050 SAEMOBJ=  -19688.3181150357
 iteration         1060 SAEMOBJ=  -19688.3139012184
 iteration         1070 SAEMOBJ=  -19688.2200016760
 iteration         1080 SAEMOBJ=  -19688.1968597214
 iteration         1090 SAEMOBJ=  -19688.2712647464
 iteration         1100 SAEMOBJ=  -19688.3720023702
 iteration         1110 SAEMOBJ=  -19688.4924690598
 iteration         1120 SAEMOBJ=  -19688.4693084105
 iteration         1130 SAEMOBJ=  -19688.4350437788
 iteration         1140 SAEMOBJ=  -19688.4897724820
 iteration         1150 SAEMOBJ=  -19688.5261776539
 iteration         1160 SAEMOBJ=  -19688.4423936266
 iteration         1170 SAEMOBJ=  -19688.5024855648
 iteration         1180 SAEMOBJ=  -19688.4921775635
 iteration         1190 SAEMOBJ=  -19688.4941076617
 iteration         1200 SAEMOBJ=  -19688.5194556680
 iteration         1210 SAEMOBJ=  -19688.5251581982
 iteration         1220 SAEMOBJ=  -19688.5197590879
 iteration         1230 SAEMOBJ=  -19688.4806140649
 iteration         1240 SAEMOBJ=  -19688.5420225876
 iteration         1250 SAEMOBJ=  -19688.5475895000
 iteration         1260 SAEMOBJ=  -19688.5289929995
 iteration         1270 SAEMOBJ=  -19688.5093730039
 iteration         1280 SAEMOBJ=  -19688.5396565915
 iteration         1290 SAEMOBJ=  -19688.4914688333
 iteration         1300 SAEMOBJ=  -19688.3834093501
 iteration         1310 SAEMOBJ=  -19688.3008599055
 iteration         1320 SAEMOBJ=  -19688.3004276725
 iteration         1330 SAEMOBJ=  -19688.2815867849
 iteration         1340 SAEMOBJ=  -19688.2740573988
 iteration         1350 SAEMOBJ=  -19688.3241011015
 iteration         1360 SAEMOBJ=  -19688.2298452081
 iteration         1370 SAEMOBJ=  -19688.3263429818
 iteration         1380 SAEMOBJ=  -19688.3554948982
 iteration         1390 SAEMOBJ=  -19688.3874986016
 iteration         1400 SAEMOBJ=  -19688.3477145957
 iteration         1410 SAEMOBJ=  -19688.2494467411
 iteration         1420 SAEMOBJ=  -19688.3406424311
 iteration         1430 SAEMOBJ=  -19688.3275869107
 iteration         1440 SAEMOBJ=  -19688.2372420392
 iteration         1450 SAEMOBJ=  -19688.2460843030
 iteration         1460 SAEMOBJ=  -19688.2254373165
 iteration         1470 SAEMOBJ=  -19688.1776494826
 iteration         1480 SAEMOBJ=  -19688.2545422339
 iteration         1490 SAEMOBJ=  -19688.3523206556
 iteration         1500 SAEMOBJ=  -19688.1844715366
 iteration         1510 SAEMOBJ=  -19688.0058934417
 iteration         1520 SAEMOBJ=  -19687.9513826406
 iteration         1530 SAEMOBJ=  -19688.0334101776
 iteration         1540 SAEMOBJ=  -19688.0381936201
 iteration         1550 SAEMOBJ=  -19687.9650493731
 iteration         1560 SAEMOBJ=  -19687.9513945509
 iteration         1570 SAEMOBJ=  -19687.9117934234
 iteration         1580 SAEMOBJ=  -19687.9123237516
 iteration         1590 SAEMOBJ=  -19687.7968733871
 iteration         1600 SAEMOBJ=  -19687.8007165770
 iteration         1610 SAEMOBJ=  -19687.7832537256
 iteration         1620 SAEMOBJ=  -19687.7700053795
 iteration         1630 SAEMOBJ=  -19687.7296129062
 iteration         1640 SAEMOBJ=  -19687.7408536442
 iteration         1650 SAEMOBJ=  -19687.7276260127
 iteration         1660 SAEMOBJ=  -19687.7198309366
 iteration         1670 SAEMOBJ=  -19687.7468580429
 iteration         1680 SAEMOBJ=  -19687.7899436601
 iteration         1690 SAEMOBJ=  -19687.8499146979
 iteration         1700 SAEMOBJ=  -19687.8146163800
 iteration         1710 SAEMOBJ=  -19687.7573997118
 iteration         1720 SAEMOBJ=  -19687.6708719811
 iteration         1730 SAEMOBJ=  -19687.6639401151
 iteration         1740 SAEMOBJ=  -19687.6208806474
 iteration         1750 SAEMOBJ=  -19687.4891923449
 iteration         1760 SAEMOBJ=  -19687.4612136719
 iteration         1770 SAEMOBJ=  -19687.4442714723
 iteration         1780 SAEMOBJ=  -19687.4679644472
 iteration         1790 SAEMOBJ=  -19687.4914500116
 iteration         1800 SAEMOBJ=  -19687.4712687674
 iteration         1810 SAEMOBJ=  -19687.5532304209
 iteration         1820 SAEMOBJ=  -19687.5894964726
 iteration         1830 SAEMOBJ=  -19687.5292219358
 iteration         1840 SAEMOBJ=  -19687.5592492452
 iteration         1850 SAEMOBJ=  -19687.5762786909
 iteration         1860 SAEMOBJ=  -19687.5098893862
 iteration         1870 SAEMOBJ=  -19687.6203103404
 iteration         1880 SAEMOBJ=  -19687.7269156143
 iteration         1890 SAEMOBJ=  -19687.6811572920
 iteration         1900 SAEMOBJ=  -19687.7121921442
 iteration         1910 SAEMOBJ=  -19687.7535043869
 iteration         1920 SAEMOBJ=  -19687.7279073482
 iteration         1930 SAEMOBJ=  -19687.7070756956
 iteration         1940 SAEMOBJ=  -19687.7577447787
 iteration         1950 SAEMOBJ=  -19687.7698371646
 iteration         1960 SAEMOBJ=  -19687.7989787757
 iteration         1970 SAEMOBJ=  -19687.8355868108
 iteration         1980 SAEMOBJ=  -19687.7990496404
 iteration         1990 SAEMOBJ=  -19687.7900491282
 iteration         2000 SAEMOBJ=  -19687.8231967026
 
 #TERM:
 STOCHASTIC PORTION WAS COMPLETED
 REDUCED STOCHASTIC PORTION WAS COMPLETED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         4.0942E-07 -8.1929E-07  5.4579E-08  2.4478E-08
 SE:             4.7450E-03  2.8819E-03  2.8503E-03  3.6986E-03
 N:                     400         400         400         400
 
 P VAL.:         9.9993E-01  9.9977E-01  9.9998E-01  9.9999E-01
 
 ETASHRINKSD(%)  7.0122E+00  3.4102E+01  4.1899E+01  2.5220E+01
 ETASHRINKVR(%)  1.3533E+01  5.6574E+01  6.6243E+01  4.4079E+01
 EBVSHRINKSD(%)  7.0112E+00  3.4102E+01  4.1898E+01  2.5222E+01
 EBVSHRINKVR(%)  1.3531E+01  5.6574E+01  6.6242E+01  4.4082E+01
 EPSSHRINKSD(%)  2.6247E+01
 EPSSHRINKVR(%)  4.5605E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2000
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    3675.75413281869     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -19687.8231967026     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -16012.0690638839     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1600
 NIND*NETA*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    2940.60330625495     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -19687.8231967026     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -16747.2198904477     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 #TERE:
 Elapsed estimation  time in seconds:   197.01
 Elapsed covariance  time in seconds:     0.43
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************           STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION (NO PRIOR)         ********************
 #OBJT:**************                        FINAL VALUE OF LIKELIHOOD FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -19687.823       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************           STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION (NO PRIOR)         ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.30E+00  3.25E+00 -6.12E-01 -2.08E-01  7.32E-01  1.14E+00  3.36E-01  1.91E-01  6.91E-01  2.30E+00  1.00E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.04E-02
 
 ETA2
+        2.20E-04  7.67E-03
 
 ETA3
+        1.46E-03 -2.16E-04  9.65E-03
 
 ETA4
+       -4.10E-04  5.31E-04  1.90E-03  9.81E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.02E-01
 
 ETA2
+        2.46E-02  8.76E-02
 
 ETA3
+        1.45E-01 -2.51E-02  9.82E-02
 
 ETA4
+       -4.05E-02  6.12E-02  1.95E-01  9.90E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************           STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION (NO PRIOR)         ********************
 ********************                          STANDARD ERROR OF ESTIMATE (S)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.82E-02  2.89E-02  1.11E-02  8.44E-03  4.80E-02  4.13E-02  1.33E-02  1.19E-02  1.05E-02  9.01E-03  3.02E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.03E-03
 
 ETA2
+        8.69E-04  1.49E-03
 
 ETA3
+        1.31E-03  1.43E-03  2.78E-03
 
 ETA4
+        1.14E-03  1.16E-03  1.96E-03  2.01E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        5.02E-03
 
 ETA2
+        9.58E-02  8.50E-03
 
 ETA3
+        1.22E-01  1.68E-01  1.42E-02
 
 ETA4
+        1.16E-01  1.30E-01  1.64E-01  1.02E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************           STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION (NO PRIOR)         ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (S)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.46E-03
 
 TH 2
+        7.67E-07  8.33E-04
 
 TH 3
+       -4.16E-04  2.15E-06  1.24E-04
 
 TH 4
+        6.18E-08 -2.33E-04 -2.95E-07  7.12E-05
 
 TH 5
+        6.33E-04  4.69E-05 -1.65E-04 -3.29E-06  2.31E-03
 
 TH 6
+       -1.67E-05  2.26E-04  6.12E-06 -5.52E-05  2.25E-05  1.70E-03
 
 TH 7
+       -1.66E-04 -1.41E-05  4.38E-05  1.33E-06 -6.27E-04 -9.24E-06  1.78E-04
 
 TH 8
+        4.95E-06 -5.64E-05 -1.38E-06  1.56E-05  3.13E-06 -4.75E-04 -2.36E-08  1.41E-04
 
 TH 9
+        1.60E-05  5.93E-05  6.06E-07 -8.68E-06  1.12E-04  5.42E-05 -3.13E-05 -7.68E-06  1.11E-04
 
 TH10
+        1.17E-05  3.58E-05  8.61E-07 -4.33E-06  1.03E-04  4.95E-05 -2.84E-05 -5.61E-06  6.66E-05  8.12E-05
 
 TH11
+       -6.84E-06 -3.00E-06  2.41E-06  1.21E-06 -6.55E-06  1.96E-06  1.98E-06  2.88E-08  3.21E-06  2.94E-06  9.13E-06
 
 OM11
+        4.01E-06 -4.99E-07 -1.31E-06 -1.16E-07  9.92E-07  1.32E-07 -3.87E-07 -7.88E-08 -6.90E-07 -1.48E-07 -4.17E-07  1.05E-06
 
 OM12
+        1.31E-06  5.88E-09 -6.16E-07 -5.59E-08  1.37E-06  4.15E-07 -2.76E-07 -1.68E-09 -4.93E-07  3.91E-08 -5.58E-07  3.79E-07
          7.55E-07
 
 OM13
+        3.36E-07  1.10E-06 -4.19E-07 -4.44E-07  1.64E-06 -6.70E-07 -4.58E-07  1.49E-07  1.01E-07  1.20E-07 -4.92E-07  5.67E-07
          3.43E-07  1.72E-06
 
 OM14
+        1.77E-06 -1.18E-06 -7.45E-07  4.34E-07  5.62E-06 -8.54E-07 -1.39E-06  2.34E-07  1.94E-07  4.27E-07 -3.17E-07  4.65E-07
          3.11E-07  1.08E-06  1.31E-06
 
 OM22
+        3.22E-06  1.33E-07 -9.54E-07 -4.78E-08  4.23E-06 -7.60E-06 -5.66E-07  2.06E-06 -8.99E-07 -7.91E-07 -1.16E-06  8.68E-08
          4.62E-07  7.32E-08  1.63E-07  2.21E-06
 
 OM23
+        1.04E-07  4.47E-08 -2.35E-07 -2.16E-07 -8.34E-07  6.07E-08  2.59E-07  1.15E-07 -1.28E-06 -9.22E-07 -7.12E-07  3.10E-07
          4.41E-07  7.91E-07  5.17E-07  5.47E-07  2.03E-06
 
 OM24
+        2.39E-06 -8.82E-07 -6.56E-07  1.48E-07  4.72E-06 -2.81E-06 -1.21E-06  1.01E-06 -4.14E-07  4.76E-07 -3.76E-07  2.10E-07
          3.76E-07  4.80E-07  5.54E-07  5.41E-07  1.11E-06  1.35E-06
 
 OM33
+        3.10E-06  9.08E-07 -1.07E-06 -4.42E-07  8.15E-06 -6.63E-06 -2.52E-06  1.23E-06 -2.79E-06 -2.54E-06 -3.17E-06  4.69E-07
          4.55E-07  1.32E-06  6.64E-07  2.70E-07  1.56E-06  7.09E-07  7.74E-06
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        2.09E-06  2.08E-06 -7.80E-07 -7.07E-07  6.42E-06  1.07E-06 -2.09E-06 -4.77E-07 -1.54E-06 -7.40E-07 -2.28E-06  4.37E-07
          3.97E-07  9.29E-07  8.74E-07  2.09E-07  1.21E-06  8.02E-07  4.37E-06  3.84E-06
 
 OM44
+       -8.10E-07  2.91E-06  8.10E-08 -6.94E-07  7.63E-06  4.88E-06 -2.37E-06 -9.94E-07  1.86E-07  2.35E-06 -1.82E-06  4.32E-07
          3.85E-07  7.65E-07  1.02E-06  9.02E-08  8.14E-07  9.69E-07  2.30E-06  3.01E-06  4.05E-06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************           STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION (NO PRIOR)         ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (S)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        3.82E-02
 
 TH 2
+        6.96E-04  2.89E-02
 
 TH 3
+       -9.77E-01  6.69E-03  1.11E-02
 
 TH 4
+        1.92E-04 -9.55E-01 -3.13E-03  8.44E-03
 
 TH 5
+        3.45E-01  3.38E-02 -3.09E-01 -8.13E-03  4.80E-02
 
 TH 6
+       -1.06E-02  1.90E-01  1.33E-02 -1.59E-01  1.14E-02  4.13E-02
 
 TH 7
+       -3.27E-01 -3.67E-02  2.95E-01  1.19E-02 -9.80E-01 -1.68E-02  1.33E-02
 
 TH 8
+        1.09E-02 -1.64E-01 -1.04E-02  1.56E-01  5.48E-03 -9.68E-01 -1.49E-04  1.19E-02
 
 TH 9
+        3.98E-02  1.95E-01  5.16E-03 -9.76E-02  2.20E-01  1.25E-01 -2.23E-01 -6.13E-02  1.05E-02
 
 TH10
+        3.39E-02  1.38E-01  8.57E-03 -5.69E-02  2.39E-01  1.33E-01 -2.37E-01 -5.23E-02  7.01E-01  9.01E-03
 
 TH11
+       -5.93E-02 -3.44E-02  7.16E-02  4.73E-02 -4.52E-02  1.57E-02  4.91E-02  8.02E-04  1.01E-01  1.08E-01  3.02E-03
 
 OM11
+        1.03E-01 -1.69E-02 -1.15E-01 -1.34E-02  2.02E-02  3.12E-03 -2.83E-02 -6.46E-03 -6.39E-02 -1.60E-02 -1.35E-01  1.03E-03
 
 OM12
+        3.94E-02  2.34E-04 -6.36E-02 -7.62E-03  3.29E-02  1.16E-02 -2.38E-02 -1.63E-04 -5.38E-02  5.00E-03 -2.13E-01  4.25E-01
          8.69E-04
 
 OM13
+        6.72E-03  2.92E-02 -2.87E-02 -4.01E-02  2.61E-02 -1.24E-02 -2.62E-02  9.55E-03  7.33E-03  1.02E-02 -1.24E-01  4.23E-01
          3.02E-01  1.31E-03
 
 OM14
+        4.07E-02 -3.58E-02 -5.85E-02  4.50E-02  1.02E-01 -1.81E-02 -9.10E-02  1.73E-02  1.61E-02  4.14E-02 -9.18E-02  3.97E-01
          3.13E-01  7.18E-01  1.14E-03
 
 OM22
+        5.66E-02  3.09E-03 -5.75E-02 -3.81E-03  5.92E-02 -1.24E-01 -2.85E-02  1.17E-01 -5.73E-02 -5.90E-02 -2.59E-01  5.69E-02
          3.58E-01  3.75E-02  9.56E-02  1.49E-03
 
 OM23
+        1.91E-03  1.09E-03 -1.48E-02 -1.79E-02 -1.22E-02  1.03E-03  1.36E-02  6.76E-03 -8.49E-02 -7.17E-02 -1.65E-01  2.12E-01
          3.56E-01  4.23E-01  3.17E-01  2.58E-01  1.43E-03
 
 OM24
+        5.39E-02 -2.63E-02 -5.06E-02  1.51E-02  8.46E-02 -5.85E-02 -7.80E-02  7.27E-02 -3.38E-02  4.55E-02 -1.07E-01  1.77E-01
          3.72E-01  3.15E-01  4.17E-01  3.12E-01  6.71E-01  1.16E-03
 
 OM33
+        2.92E-02  1.13E-02 -3.44E-02 -1.88E-02  6.10E-02 -5.77E-02 -6.80E-02  3.73E-02 -9.52E-02 -1.01E-01 -3.77E-01  1.64E-01
          1.88E-01  3.62E-01  2.09E-01  6.52E-02  3.92E-01  2.19E-01  2.78E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        2.80E-02  3.67E-02 -3.57E-02 -4.27E-02  6.82E-02  1.32E-02 -8.01E-02 -2.05E-02 -7.47E-02 -4.19E-02 -3.84E-01  2.18E-01
          2.33E-01  3.62E-01  3.90E-01  7.17E-02  4.32E-01  3.52E-01  8.02E-01  1.96E-03
 
 OM44
+       -1.05E-02  5.02E-02  3.61E-03 -4.08E-02  7.90E-02  5.87E-02 -8.83E-02 -4.16E-02  8.77E-03  1.29E-01 -2.99E-01  2.10E-01
          2.20E-01  2.90E-01  4.41E-01  3.01E-02  2.83E-01  4.14E-01  4.10E-01  7.62E-01  2.01E-03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************           STOCHASTIC APPROXIMATION EXPECTATION-MAXIMIZATION (NO PRIOR)         ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (S)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.70E+04
 
 TH 2
+        0.00E+00  1.71E+04
 
 TH 3
+        5.60E+04  0.00E+00  1.94E+05
 
 TH 4
+        0.00E+00  5.53E+04  0.00E+00  1.94E+05
 
 TH 5
+       -2.05E+03  0.00E+00 -5.68E+03  0.00E+00  1.17E+04
 
 TH 6
+        0.00E+00 -3.62E+03  0.00E+00 -1.16E+04  0.00E+00  1.16E+04
 
 TH 7
+       -5.68E+03  0.00E+00 -1.72E+04  0.00E+00  4.05E+04  0.00E+00  1.48E+05
 
 TH 8
+        0.00E+00 -1.16E+04  0.00E+00 -3.87E+04  0.00E+00  3.85E+04  0.00E+00  1.35E+05
 
 TH 9
+       -1.64E+03 -3.68E+03 -5.65E+03 -1.06E+04  5.14E+02  3.50E+00  2.69E+03  1.84E+02  1.95E+04
 
 TH10
+       -1.24E+03 -1.39E+02 -4.66E+03 -8.39E+02 -6.93E+02 -3.21E+03 -9.60E+02 -1.05E+04 -1.47E+04  2.78E+04
 
 TH11
+       -7.50E+02 -4.75E+02 -4.33E+03 -2.89E+03 -1.22E+03 -7.51E+02 -5.22E+03 -3.25E+03  1.35E+02 -4.42E+03  1.48E+05
 
 OM11
+       -9.84E+03  1.40E+04 -2.07E+04  4.71E+04  9.70E+03 -5.68E+02  3.22E+04 -2.24E+03  7.13E+03 -4.08E+03  1.60E+04  1.40E+06
 
 OM12
+        2.01E+04 -1.22E+03  6.63E+04 -5.03E+03 -4.47E+03 -9.47E+03 -1.61E+04 -2.54E+04  5.76E+03 -1.00E+04  3.71E+04 -5.82E+05
          2.08E+06
 
 OM13
+        6.51E+03  9.67E+03  1.99E+04  4.81E+04  5.67E+03 -1.76E+03  1.52E+04 -7.04E+03 -6.90E+03 -4.88E+03  3.80E+03 -2.44E+05
         -2.99E+04  1.65E+06
 
 OM14
+        9.74E+03 -1.78E+04  3.24E+04 -8.06E+04 -1.83E+04  9.92E+03 -5.49E+04  3.21E+04 -5.42E+03  1.13E+04 -3.53E+04 -2.06E+05
         -6.42E+04 -1.29E+06  2.26E+06
 
 OM22
+       -1.41E+03 -5.04E+03 -4.67E+03 -1.35E+04 -1.16E+04  5.28E+03 -3.99E+04  9.95E+03 -3.23E+02  2.62E+03  7.58E+04  6.72E+04
         -3.12E+05  8.93E+04 -4.26E+04  6.28E+05
 
 OM23
+        7.62E+03 -7.18E+00  2.21E+04 -1.21E+03 -6.39E+02 -8.66E+03 -7.48E+03 -2.22E+04 -1.55E+03  1.00E+04 -3.04E+03 -1.92E+04
         -9.86E+04 -3.95E+05  3.07E+05 -5.33E+04  1.22E+06
 
 OM24
+       -1.76E+04  1.12E+04 -5.31E+04  3.59E+04  4.19E+03  4.36E+02  2.00E+04 -8.69E+03  7.76E+03 -1.18E+04 -4.85E+04  1.31E+05
         -2.30E+05  2.40E+05 -4.79E+05 -1.81E+05 -9.18E+05  1.86E+06
 
 OM33
+       -4.16E+02 -6.72E+03 -3.92E+03 -2.65E+04 -4.62E+03  5.76E+03 -1.26E+04  1.37E+04  3.72E+03 -4.52E+03  4.28E+04  1.69E+04
         -3.40E+04 -2.88E+05  3.06E+05  3.35E+04  1.61E+04 -4.07E+04  5.76E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+       -7.39E+03  8.53E+03 -1.74E+04  3.57E+04  1.05E+04 -1.26E+03  3.61E+04  2.30E+03 -8.51E+03  2.65E+04  2.29E+02 -6.89E+02
          4.34E+04  3.43E+05 -4.63E+05 -3.24E+04 -3.58E+05  2.92E+05 -8.95E+05  2.22E+06
 
 OM44
+        7.18E+03 -5.09E+03  1.65E+04 -1.54E+04 -1.64E+03 -4.37E+03 -5.59E+03 -1.40E+04  1.30E+04 -3.24E+04  5.69E+04 -2.77E+04
         -1.66E+04 -2.46E+04 -9.29E+04  8.96E+04  2.37E+05 -3.86E+05  3.42E+05 -1.10E+06  9.94E+05
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 
 
 #TBLN:      4
 #METH: Objective Function Evaluation by Importance Sampling (No Prior)
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            2208
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      8
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     8
 NOPRIOR SETTING (NOPRIOR):                 ON
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example2.ext
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
 GRADIENT/GIBBS PATTERN (GRD):              DDDDDDDDDDS
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  3
 CONVERGENCE INTERVAL (CINTERVAL):          1
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 ITERATIONS (NITER):                        5
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       123334
 MC SAMPLES PER SUBJECT (ISAMPLE):          3000
 RANDOM SAMPLING METHOD (RANMETHOD):        3U
 EXPECTATION ONLY (EONLY):                  1
 PROPOSAL DENSITY SCALING RANGE
              (ISCALE_MIN, ISCALE_MAX):     0.100000000000000       ,10.0000000000000
 SAMPLE ACCEPTANCE RATE (IACCEPT):          0.400000000000000
 LONG TAIL SAMPLE ACCEPT. RATE (IACCEPTL):   0.00000000000000
 T-DIST. PROPOSAL DENSITY (DF):             0
 NO. ITERATIONS FOR MAP (MAPITER):          0
 INTERVAL ITER. FOR MAP (MAPINTER):         0
 MAP COVARIANCE/MODE SETTING (MAPCOV):      1
 Gradient Quick Value (GRDQ):               0.00000000000000

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8   9  10
 THETAS THAT ARE SIGMA-LIKE:
  11
 
 MONITORING OF SEARCH:

 iteration            0 OBJ=  -10782.2889938195 eff.=    2975. Smpl.=    3000. Fit.= 0.98016
 iteration            1 OBJ=  -10782.3764633088 eff.=    1201. Smpl.=    3000. Fit.= 0.90045
 iteration            2 OBJ=  -10782.7883020672 eff.=    1195. Smpl.=    3000. Fit.= 0.89986
 iteration            3 OBJ=  -10781.9935886524 eff.=    1199. Smpl.=    3000. Fit.= 0.90036
 iteration            4 OBJ=  -10781.8220380917 eff.=    1196. Smpl.=    3000. Fit.= 0.90030
 iteration            5 OBJ=  -10783.0603980877 eff.=    1198. Smpl.=    3000. Fit.= 0.90016
 
 #TERM:
 EXPECTATION ONLY PROCESS COMPLETED


 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -9.0969E-05 -1.2242E-04 -1.6044E-04 -2.1097E-04
 SE:             4.7407E-03  2.8957E-03  2.8548E-03  3.6929E-03
 N:                     400         400         400         400
 
 P VAL.:         9.8469E-01  9.6628E-01  9.5518E-01  9.5444E-01
 
 ETASHRINKSD(%)  7.0955E+00  3.3786E+01  4.1807E+01  2.5335E+01
 ETASHRINKVR(%)  1.3688E+01  5.6157E+01  6.6135E+01  4.4251E+01
 EBVSHRINKSD(%)  6.9839E+00  3.4234E+01  4.2026E+01  2.5161E+01
 EBVSHRINKVR(%)  1.3480E+01  5.6748E+01  6.6391E+01  4.3991E+01
 EPSSHRINKSD(%)  2.6291E+01
 EPSSHRINKVR(%)  4.5670E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2000
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    3675.75413281869     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -10783.0603980877     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -7107.30626526901     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1600
  
 #TERE:
 Elapsed estimation  time in seconds:    27.76
 Elapsed covariance  time in seconds:    47.98
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************         OBJECTIVE FUNCTION EVALUATION BY IMPORTANCE SAMPLING (NO PRIOR)        ********************
 #OBJT:**************                        FINAL VALUE OF OBJECTIVE FUNCTION                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -10783.060       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************         OBJECTIVE FUNCTION EVALUATION BY IMPORTANCE SAMPLING (NO PRIOR)        ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.30E+00  3.25E+00 -6.12E-01 -2.08E-01  7.32E-01  1.14E+00  3.36E-01  1.91E-01  6.91E-01  2.30E+00  1.00E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.04E-02
 
 ETA2
+        2.20E-04  7.67E-03
 
 ETA3
+        1.46E-03 -2.16E-04  9.65E-03
 
 ETA4
+       -4.10E-04  5.31E-04  1.90E-03  9.81E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.02E-01
 
 ETA2
+        2.46E-02  8.76E-02
 
 ETA3
+        1.45E-01 -2.51E-02  9.82E-02
 
 ETA4
+       -4.05E-02  6.12E-02  1.95E-01  9.90E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************         OBJECTIVE FUNCTION EVALUATION BY IMPORTANCE SAMPLING (NO PRIOR)        ********************
 ********************                          STANDARD ERROR OF ESTIMATE (R)                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.28E-02  2.87E-02  9.56E-03  8.32E-03  3.90E-02  3.57E-02  1.13E-02  1.03E-02  1.04E-02  8.59E-03  2.83E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        9.87E-04
 
 ETA2
+        8.27E-04  1.34E-03
 
 ETA3
+        1.30E-03  1.43E-03  3.07E-03
 
 ETA4
+        1.03E-03  1.10E-03  2.20E-03  1.97E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        4.83E-03
 
 ETA2
+        9.13E-02  7.66E-03
 
 ETA3
+        1.16E-01  1.68E-01  1.56E-02
 
 ETA4
+        1.04E-01  1.24E-01  1.83E-01  9.95E-03
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************         OBJECTIVE FUNCTION EVALUATION BY IMPORTANCE SAMPLING (NO PRIOR)        ********************
 ********************                        COVARIANCE MATRIX OF ESTIMATE (R)                       ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.07E-03
 
 TH 2
+        1.58E-05  8.22E-04
 
 TH 3
+       -3.03E-04 -2.86E-06  9.13E-05
 
 TH 4
+       -1.01E-06 -2.29E-04  8.73E-08  6.92E-05
 
 TH 5
+        2.90E-04  2.54E-05 -7.83E-05 -2.55E-06  1.52E-03
 
 TH 6
+        2.82E-05  2.91E-04 -5.45E-06 -7.58E-05  3.96E-05  1.27E-03
 
 TH 7
+       -8.09E-05 -8.08E-06  2.28E-05  8.38E-07 -4.26E-04 -1.24E-05  1.27E-04
 
 TH 8
+       -4.66E-06 -7.78E-05  9.12E-07  2.23E-05 -5.58E-06 -3.54E-04  1.84E-06  1.07E-04
 
 TH 9
+        4.43E-05  4.56E-05 -9.32E-06 -4.70E-06  5.83E-05  7.00E-05 -1.85E-05 -1.20E-05  1.09E-04
 
 TH10
+        2.76E-05  3.19E-05 -5.05E-06 -3.61E-06  5.25E-05  5.92E-05 -1.54E-05 -9.36E-06  6.26E-05  7.37E-05
 
 TH11
+        1.61E-06 -2.77E-07 -4.50E-07  1.48E-07  4.36E-06 -2.27E-06 -1.26E-06  7.99E-07  1.29E-06  9.63E-07  7.99E-06
 
 OM11
+        2.05E-07 -2.12E-07 -9.85E-09  9.33E-08  1.31E-07  3.33E-07  4.10E-08 -2.96E-08  2.57E-07  2.20E-07 -3.93E-07  9.74E-07
 
 OM12
+        8.53E-08 -2.27E-07  2.18E-08  1.19E-07 -7.26E-08 -2.78E-07  1.49E-07  1.67E-07  2.13E-07  1.79E-07 -4.28E-07  2.83E-07
          6.84E-07
 
 OM13
+        3.16E-07 -2.37E-08  5.58E-08  4.15E-08  5.24E-07  1.02E-06 -3.53E-08 -1.59E-07  5.62E-07  5.51E-07 -7.35E-07  5.77E-07
          2.28E-07  1.68E-06
 
 OM14
+        8.64E-07 -1.50E-06 -1.96E-07  5.26E-07  1.64E-07  7.94E-07  1.40E-08 -1.36E-07  5.42E-07  3.33E-07 -6.08E-07  3.95E-07
          2.25E-07  9.42E-07  1.06E-06
 
 OM22
+       -4.99E-07 -2.32E-07  1.77E-07  8.49E-08 -1.63E-06 -1.17E-06  6.46E-07  4.17E-07 -4.27E-07 -3.20E-07 -1.22E-06  8.38E-08
          3.63E-07  6.12E-08  7.31E-08  1.80E-06
 
 OM23
+        5.49E-07  2.53E-07 -7.39E-08 -9.62E-10  1.84E-06 -3.66E-07 -2.50E-07  2.32E-07  1.10E-06  6.89E-07 -6.53E-07  1.96E-07
          3.69E-07  7.34E-07  4.07E-07  3.53E-07  2.04E-06
 
 OM24
+        2.63E-07 -2.09E-07 -2.96E-08  1.04E-07  1.44E-06 -2.33E-06 -3.09E-07  7.88E-07  5.94E-07  4.10E-07 -4.36E-07  1.10E-07
          2.51E-07  3.27E-07  3.55E-07  4.07E-07  1.00E-06  1.21E-06
 
 OM33
+       -3.19E-07  1.67E-06  2.81E-07 -2.93E-07 -4.44E-07  5.04E-06  2.24E-07 -1.14E-06  7.54E-07  8.06E-07 -3.54E-06  4.90E-07
          3.47E-07  2.10E-06  1.17E-06  3.07E-07  1.80E-06  7.44E-07  9.44E-06
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+       -2.98E-07  1.98E-06  2.09E-07 -4.23E-07 -9.57E-07  5.43E-06  3.07E-07 -1.35E-06  9.32E-07  6.28E-07 -2.54E-06  3.31E-07
          2.57E-07  1.35E-06  1.02E-06  1.78E-07  1.32E-06  7.09E-07  5.81E-06  4.82E-06
 
 OM44
+       -3.91E-07  1.89E-06  1.91E-07 -4.32E-07 -1.35E-06  5.41E-06  3.82E-07 -1.42E-06  6.88E-07  5.72E-07 -1.94E-06  2.20E-07
          1.84E-07  8.01E-07  8.33E-07  1.23E-07  7.83E-07  7.19E-07  3.36E-06  3.50E-06  3.89E-06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************         OBJECTIVE FUNCTION EVALUATION BY IMPORTANCE SAMPLING (NO PRIOR)        ********************
 ********************                        CORRELATION MATRIX OF ESTIMATE (R)                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        3.28E-02
 
 TH 2
+        1.68E-02  2.87E-02
 
 TH 3
+       -9.68E-01 -1.04E-02  9.56E-03
 
 TH 4
+       -3.69E-03 -9.59E-01  1.10E-03  8.32E-03
 
 TH 5
+        2.27E-01  2.27E-02 -2.10E-01 -7.84E-03  3.90E-02
 
 TH 6
+        2.41E-02  2.84E-01 -1.60E-02 -2.55E-01  2.85E-02  3.57E-02
 
 TH 7
+       -2.19E-01 -2.50E-02  2.12E-01  8.93E-03 -9.69E-01 -3.08E-02  1.13E-02
 
 TH 8
+       -1.38E-02 -2.63E-01  9.24E-03  2.59E-01 -1.39E-02 -9.59E-01  1.58E-02  1.03E-02
 
 TH 9
+        1.30E-01  1.52E-01 -9.35E-02 -5.41E-02  1.43E-01  1.88E-01 -1.57E-01 -1.12E-01  1.04E-02
 
 TH10
+        9.81E-02  1.30E-01 -6.16E-02 -5.05E-02  1.57E-01  1.93E-01 -1.59E-01 -1.06E-01  6.99E-01  8.59E-03
 
 TH11
+        1.74E-02 -3.42E-03 -1.67E-02  6.29E-03  3.95E-02 -2.25E-02 -3.94E-02  2.74E-02  4.37E-02  3.97E-02  2.83E-03
 
 OM11
+        6.33E-03 -7.47E-03 -1.04E-03  1.14E-02  3.41E-03  9.45E-03  3.68E-03 -2.90E-03  2.50E-02  2.60E-02 -1.41E-01  9.87E-04
 
 OM12
+        3.15E-03 -9.56E-03  2.76E-03  1.73E-02 -2.25E-03 -9.42E-03  1.59E-02  1.95E-02  2.46E-02  2.52E-02 -1.83E-01  3.47E-01
          8.27E-04
 
 OM13
+        7.43E-03 -6.36E-04  4.50E-03  3.84E-03  1.03E-02  2.19E-02 -2.42E-03 -1.19E-02  4.15E-02  4.95E-02 -2.00E-01  4.50E-01
          2.12E-01  1.30E-03
 
 OM14
+        2.57E-02 -5.09E-02 -2.00E-02  6.15E-02  4.09E-03  2.17E-02  1.20E-03 -1.28E-02  5.06E-02  3.78E-02 -2.09E-01  3.90E-01
          2.65E-01  7.06E-01  1.03E-03
 
 OM22
+       -1.14E-02 -6.04E-03  1.38E-02  7.61E-03 -3.11E-02 -2.44E-02  4.27E-02  3.01E-02 -3.05E-02 -2.78E-02 -3.21E-01  6.33E-02
          3.27E-01  3.51E-02  5.30E-02  1.34E-03
 
 OM23
+        1.17E-02  6.18E-03 -5.41E-03 -8.10E-05  3.30E-02 -7.19E-03 -1.56E-02  1.58E-02  7.40E-02  5.62E-02 -1.62E-01  1.39E-01
          3.12E-01  3.96E-01  2.77E-01  1.84E-01  1.43E-03
 
 OM24
+        7.28E-03 -6.62E-03 -2.81E-03  1.13E-02  3.36E-02 -5.94E-02 -2.49E-02  6.93E-02  5.17E-02  4.34E-02 -1.40E-01  1.01E-01
          2.76E-01  2.29E-01  3.14E-01  2.76E-01  6.37E-01  1.10E-03
 
 OM33
+       -3.17E-03  1.89E-02  9.56E-03 -1.14E-02 -3.70E-03  4.59E-02  6.47E-03 -3.59E-02  2.35E-02  3.06E-02 -4.08E-01  1.62E-01
          1.36E-01  5.27E-01  3.71E-01  7.46E-02  4.11E-01  2.20E-01  3.07E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+       -4.15E-03  3.15E-02  9.97E-03 -2.32E-02 -1.12E-02  6.93E-02  1.24E-02 -5.97E-02  4.07E-02  3.33E-02 -4.09E-01  1.53E-01
          1.41E-01  4.73E-01  4.52E-01  6.06E-02  4.20E-01  2.93E-01  8.62E-01  2.20E-03
 
 OM44
+       -6.05E-03  3.34E-02  1.01E-02 -2.63E-02 -1.75E-02  7.69E-02  1.72E-02 -6.96E-02  3.35E-02  3.38E-02 -3.48E-01  1.13E-01
          1.13E-01  3.13E-01  4.11E-01  4.66E-02  2.78E-01  3.31E-01  5.54E-01  8.09E-01  1.97E-03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************         OBJECTIVE FUNCTION EVALUATION BY IMPORTANCE SAMPLING (NO PRIOR)        ********************
 ********************                    INVERSE COVARIANCE MATRIX OF ESTIMATE (R)                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.57E+04
 
 TH 2
+        0.00E+00  1.84E+04
 
 TH 3
+        5.18E+04  0.00E+00  1.82E+05
 
 TH 4
+        0.00E+00  6.05E+04  0.00E+00  2.14E+05
 
 TH 5
+       -2.19E+03  0.00E+00 -6.97E+03  0.00E+00  1.14E+04
 
 TH 6
+        0.00E+00 -3.48E+03  0.00E+00 -1.12E+04  0.00E+00  1.18E+04
 
 TH 7
+       -6.97E+03  0.00E+00 -2.41E+04  0.00E+00  3.81E+04  0.00E+00  1.36E+05
 
 TH 8
+        0.00E+00 -1.12E+04  0.00E+00 -3.92E+04  0.00E+00  3.86E+04  0.00E+00  1.37E+05
 
 TH 9
+       -1.36E+03 -3.94E+03 -3.65E+03 -1.27E+04  1.14E+03 -5.46E+02  4.41E+03 -1.16E+03  1.95E+04
 
 TH10
+       -1.05E+03 -3.10E+02 -3.80E+03 -9.52E+02 -7.20E+02 -3.14E+03 -1.37E+03 -9.65E+03 -1.48E+04  2.79E+04
 
 TH11
+       -2.74E+02 -2.94E+02 -7.88E+02 -1.17E+03 -9.28E+02 -1.08E+03 -2.67E+03 -4.16E+03 -1.19E+03 -9.42E+02  1.76E+05
 
 OM11
+        6.91E+02 -1.91E+02  2.42E+03 -5.64E+02 -1.28E+03 -4.21E+02 -4.31E+03 -1.19E+03 -1.22E+03  2.59E+02  3.37E+04  1.45E+06
 
 OM12
+       -1.56E+03 -1.81E+03 -5.95E+03 -5.67E+03 -2.94E+03 -1.98E+03 -1.11E+04 -7.83E+03  9.99E+02 -2.49E+03  1.91E+04 -4.93E+05
          2.03E+06
 
 OM13
+       -5.63E+03  4.55E+03 -2.24E+04  2.43E+04  4.35E+02 -1.27E+03  2.00E+03 -6.76E+03  4.03E+03 -5.87E+03 -3.07E+04 -4.60E+05
          1.83E+05  1.74E+06
 
 OM14
+        3.15E+03 -3.52E+03  1.67E+04 -3.09E+04 -1.29E+03  5.81E+02 -5.68E+03  6.27E+03 -6.25E+03  5.65E+03  2.58E+04 -1.16E+05
         -3.02E+05 -1.26E+06  2.38E+06
 
 OM22
+       -7.81E+02 -1.19E+03 -2.64E+03 -4.12E+03 -2.46E+03 -2.31E+03 -1.10E+04 -7.83E+03  1.64E+03  1.87E+03  1.10E+05  4.55E+04
         -3.10E+05  1.95E+03  6.93E+04  7.32E+05
 
 OM23
+        2.23E+03  5.27E+02  8.63E+03 -6.32E+02 -6.90E+03  1.46E+03 -2.37E+04  7.77E+03 -6.55E+03  1.92E+03 -1.26E+04  8.43E+04
         -2.78E+05 -3.66E+05  3.23E+05  2.58E+04  1.11E+06
 
 OM24
+       -7.66E+02  1.58E+03 -4.25E+03  9.47E+03  2.56E+03  3.84E+02  1.27E+04 -1.11E+04  1.45E+02 -2.93E+03 -2.60E+04  1.77E+04
         -1.74E+04  2.52E+05 -4.51E+05 -2.29E+05 -8.25E+05  1.71E+06
 
 OM33
+        3.29E+01 -1.54E+03  6.32E+02 -7.57E+03 -2.77E+02 -7.51E+02 -4.83E+02 -3.00E+03  3.90E+03 -3.65E+03  5.03E+04  4.12E+04
         -1.65E+04 -3.10E+05  2.42E+05  2.02E+03  1.49E+04 -1.79E+04  6.06E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+       -4.72E+02  6.66E+02 -2.78E+03  6.70E+03  2.04E+03 -6.10E+02  6.58E+03 -1.44E+03 -7.28E+03  8.09E+03  5.01E+03  1.57E+04
          6.17E+04  2.08E+05 -3.84E+05  1.99E+04 -3.16E+05  2.45E+05 -9.22E+05  2.23E+06
 
 OM44
+       -1.47E+02 -5.80E+02 -6.96E+02 -1.54E+03 -4.43E+02 -7.60E+02 -2.77E+03  1.71E+03  3.93E+03 -5.03E+03  4.19E+04  6.44E+03
         -3.23E+03  2.11E+04 -6.36E+04  4.59E+04  2.09E+05 -3.21E+05  3.42E+05 -1.16E+06  1.05E+06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 
 
 #TBLN:      5
 #METH: MCMC Bayesian Analysis
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            2208
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      8
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     8
 NOPRIOR SETTING (NOPRIOR):                 OFF
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example2.TXT
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
 GRADIENT/GIBBS PATTERN (GRD):              DDDDDDDDDDS
 AUTOMATIC SETTING FEATURE (AUTO):          OFF
 CONVERGENCE TYPE (CTYPE):                  3
 KEEP ITERATIONS (THIN):            1
 CONVERGENCE INTERVAL (CINTERVAL):          100
 CONVERGENCE ITERATIONS (CITER):            10
 CONVERGENCE ALPHA ERROR (CALPHA):          5.000000000000000E-02
 BURN-IN ITERATIONS (NBURN):                10000
 ITERATIONS (NITER):                        3000
 ANEAL SETTING (CONSTRAIN):                 1
 STARTING SEED FOR MC METHODS (SEED):       123334
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
 MASS/IMP./POST. MATRIX REFRESH SETTING (MASSREST):      -1
 METROPOLIS HASTINGS POPULATION SAMPLING FOR NON-GIBBS
 SAMPLED THETAS AND SIGMAS:
 PROPOSAL DENSITY SCALING RANGE
              (PSCALE_MIN, PSCALE_MAX):   1.000000000000000E-02   ,1000.00000000000
 SAMPLE ACCEPTANCE RATE (PACCEPT):                       0.500000000000000
 SAMPLES FOR GLOBAL SEARCH KERNEL (PSAMPLE_M1):          1
 SAMPLES FOR LOCAL SEARCH KERNEL (PSAMPLE_M2):           1
 SAMPLES FOR LOCAL UNIVARIATE KERNEL (PSAMPLE_M3):       1
 METROPOLIS HASTINGS POPULATION SAMPLING FOR NON-GIBBS
 SAMPLED OMEGAS:
 SAMPLE ACCEPTANCE RATE (OACCEPT):                       0.500000000000000
 SAMPLES FOR GLOBAL SEARCH KERNEL (OSAMPLE_M1):          -1
 SAMPLES FOR LOCAL SEARCH KERNEL (OSAMPLE_M2):           10
 SAMPLES FOR LOCAL UNIVARIATE SEARCH KERNEL (OSAMPLE_M3):10
 DEG. FR. FOR T DIST.  PRIOR FOR THETAS (TTDF):        0.00000000000000
 DEG. FR. FOR LKJ CORRELATION PRIOR FOR OMEGAS (OLKJDF): 0.00000000000000
 WEIGHT FACTOR FOR STD PRIOR FOR OMEGAS (OVARF): 1.00000000000000
 DEG. FR. FOR LKJ CORRELATION PRIOR FOR SIGMAS (SLKJDF): 0.00000000000000
 WEIGHT FACTOR FOR STD PRIOR FOR SIGMAS (SVARF): 1.00000000000000

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 EM/BAYES SETUP:
 THETAS THAT ARE MU MODELED:
   1   2   3   4   5   6   7   8   9  10
 THETAS THAT ARE GIBBS SAMPLED:
   1   2   3   4   5   6   7   8   9  10
 THETAS THAT ARE METROPOLIS-HASTINGS SAMPLED:
  11
 SIGMAS THAT ARE GIBBS SAMPLED:
 
 SIGMAS THAT ARE METROPOLIS-HASTINGS SAMPLED:
 
 OMEGAS ARE GIBBS SAMPLED
 
 MONITORING OF SEARCH:

 Burn-in Mode
 iteration       -10000 MCMCOBJ=   -19197.7358618528     
 iteration        -9900 MCMCOBJ=   -19312.1236912377     
 iteration        -9800 MCMCOBJ=   -19061.3241625541     
 iteration        -9700 MCMCOBJ=   -19402.9342422065     
 iteration        -9600 MCMCOBJ=   -19014.5294708860     
 iteration        -9500 MCMCOBJ=   -19671.6815482652     
 iteration        -9400 MCMCOBJ=   -19286.7537262467     
 iteration        -9300 MCMCOBJ=   -19212.5868526989     
 iteration        -9200 MCMCOBJ=   -19451.9132627879     
 iteration        -9100 MCMCOBJ=   -19553.7236672918     
 iteration        -9000 MCMCOBJ=   -19059.2963524948     
 Convergence achieved
 Sampling Mode
 iteration            0 MCMCOBJ=   -19050.1316872174     
 iteration          100 MCMCOBJ=   -19150.3602277807     
 iteration          200 MCMCOBJ=   -19159.4430292719     
 iteration          300 MCMCOBJ=   -19288.8409713675     
 iteration          400 MCMCOBJ=   -19679.9157581953     
 iteration          500 MCMCOBJ=   -19227.6999513478     
 iteration          600 MCMCOBJ=   -19408.4899157109     
 iteration          700 MCMCOBJ=   -19299.3694446725     
 iteration          800 MCMCOBJ=   -19326.0284863544     
 iteration          900 MCMCOBJ=   -19126.1774319331     
 iteration         1000 MCMCOBJ=   -19086.8594324392     
 iteration         1100 MCMCOBJ=   -19406.2010019827     
 iteration         1200 MCMCOBJ=   -19120.9182459036     
 iteration         1300 MCMCOBJ=   -19321.3257225389     
 iteration         1400 MCMCOBJ=   -19255.4457119944     
 iteration         1500 MCMCOBJ=   -19242.0407219761     
 iteration         1600 MCMCOBJ=   -19265.4681004257     
 iteration         1700 MCMCOBJ=   -19101.8469244654     
 iteration         1800 MCMCOBJ=   -19234.4991753259     
 iteration         1900 MCMCOBJ=   -19253.4112755652     
 iteration         2000 MCMCOBJ=   -19155.6278537384     
 iteration         2100 MCMCOBJ=   -19291.5552346246     
 iteration         2200 MCMCOBJ=   -19118.2704864110     
 iteration         2300 MCMCOBJ=   -19322.6751055317     
 iteration         2400 MCMCOBJ=   -19157.1251588791     
 iteration         2500 MCMCOBJ=   -19279.1414246077     
 iteration         2600 MCMCOBJ=   -19240.3209458496     
 iteration         2700 MCMCOBJ=   -19318.7836782477     
 iteration         2800 MCMCOBJ=   -19054.0828909250     
 iteration         2900 MCMCOBJ=   -19253.6562398997     
 iteration         3000 MCMCOBJ=   -19388.5075678697     
 BURN-IN WAS COMPLETED
 
 #TERM:
 STATISTICAL PORTION WAS COMPLETED
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2000
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    3675.75413281869     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -19223.7267504576     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -15547.9726176389     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1600
 NIND*NETA*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    2940.60330625495     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -19223.7267504576     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -16283.1234442027     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 PRIOR CONSTANT TO OBJECTIVE FUNCTION:    70.3639128125257     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -19223.7267504576     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -19153.3628376451     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 #TERE:
 Elapsed estimation  time in seconds:   191.69
 Elapsed covariance  time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 #OBJT:**************                       AVERAGE VALUE OF LIKELIHOOD FUNCTION                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -19223.727       **************************************************
 #OBJS:********************************************      131.392 (STD) **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.30E+00  3.25E+00 -6.12E-01 -2.08E-01  7.35E-01  1.14E+00  3.35E-01  1.92E-01  6.92E-01  2.30E+00 -3.72E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.03E-02
 
 ETA2
+        3.72E-06  7.79E-03
 
 ETA3
+        9.21E-04 -5.24E-04  8.49E-03
 
 ETA4
+       -8.97E-04  2.76E-04  7.77E-04  8.99E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.01E-01
 
 ETA2
+       -2.84E-03  8.79E-02
 
 ETA3
+        9.36E-02 -7.38E-02  9.15E-02
 
 ETA4
+       -9.88E-02  2.79E-02  6.69E-02  9.44E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 ********************                STANDARD ERROR OF ESTIMATE (From Sample Variance)               ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.21E-02  2.86E-02  9.36E-03  8.32E-03  3.81E-02  3.47E-02  1.10E-02  1.01E-02  1.01E-02  8.33E-03  9.51E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        9.47E-04
 
 ETA2
+        7.96E-04  1.30E-03
 
 ETA3
+        1.10E-03  1.31E-03  2.00E-03
 
 ETA4
+        9.19E-04  1.05E-03  1.47E-03  1.65E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        4.65E-03
 
 ETA2
+        8.93E-02  7.38E-03
 
 ETA3
+        1.12E-01  1.62E-01  1.09E-02
 
 ETA4
+        9.91E-02  1.26E-01  1.59E-01  8.61E-03
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        0.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 ********************               COVARIANCE MATRIX OF ESTIMATE (From Sample Variance)             ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.03E-03
 
 TH 2
+        1.33E-05  8.18E-04
 
 TH 3
+       -2.90E-04 -1.94E-06  8.75E-05
 
 TH 4
+       -1.50E-06 -2.28E-04  1.59E-07  6.92E-05
 
 TH 5
+        2.47E-04  2.24E-05 -6.46E-05 -8.72E-07  1.45E-03
 
 TH 6
+        2.67E-05  2.58E-04 -5.34E-06 -6.61E-05  3.34E-05  1.20E-03
 
 TH 7
+       -6.76E-05 -7.05E-06  1.87E-05  3.91E-07 -4.06E-04 -1.15E-05  1.21E-04
 
 TH 8
+       -5.18E-06 -7.06E-05  1.01E-06  2.01E-05 -4.44E-06 -3.37E-04  1.77E-06  1.03E-04
 
 TH 9
+        3.47E-05  3.88E-05 -7.75E-06 -2.59E-06  5.89E-05  7.67E-05 -1.73E-05 -1.50E-05  1.02E-04
 
 TH10
+        2.33E-05  2.77E-05 -4.70E-06 -2.31E-06  4.36E-05  6.38E-05 -1.16E-05 -1.17E-05  5.59E-05  6.94E-05
 
 TH11
+        5.65E-05  7.07E-05 -1.68E-05 -1.58E-05 -1.80E-04 -7.22E-05  5.15E-05  3.47E-05 -1.63E-07 -3.81E-05  9.05E-03
 
 OM11
+        2.61E-08 -8.10E-07  2.91E-08  2.60E-07  8.44E-07 -2.41E-07 -1.44E-07  7.90E-08  4.46E-07 -5.78E-08  8.91E-07  8.96E-07
 
 OM12
+        8.32E-07  3.03E-08 -2.14E-07  2.13E-08  1.67E-06  1.42E-07 -3.94E-07  1.20E-07  2.85E-07  3.26E-07  6.73E-08  2.41E-07
          6.33E-07
 
 OM13
+        4.15E-07 -1.65E-07  2.20E-08 -8.90E-09  8.09E-07 -7.52E-07 -2.45E-08  2.66E-07  4.24E-07  4.70E-09 -7.37E-06  4.22E-07
          1.75E-07  1.21E-06
 
 OM14
+        1.28E-06  6.55E-08 -2.41E-07  1.99E-09  1.14E-06  4.59E-08 -1.60E-07 -2.00E-08  5.82E-07  8.96E-08 -3.88E-06  2.95E-07
          1.62E-07  6.35E-07  8.45E-07
 
 OM22
+        3.74E-08  8.59E-07  4.97E-09 -1.12E-07  1.24E-06 -5.10E-07 -2.79E-07  3.85E-07  4.61E-07  3.43E-07  8.56E-07  1.58E-08
          2.66E-07 -4.04E-08  2.43E-08  1.69E-06
 
 OM23
+        5.44E-07 -1.92E-06 -6.63E-08  4.78E-07  6.22E-07 -1.32E-06  6.20E-08  4.98E-07  7.19E-07  9.81E-08  1.86E-06  1.33E-07
          3.06E-07  4.38E-07  1.98E-07  2.21E-07  1.71E-06
 
 OM24
+        6.85E-07 -6.83E-07 -7.95E-08  1.70E-07  1.97E-06 -6.43E-07 -5.29E-07  2.75E-07  8.41E-07  3.01E-07  2.42E-06  9.04E-08
          1.52E-07  1.99E-07  2.53E-07  2.93E-07  7.93E-07  1.11E-06
 
 OM33
+        2.21E-06  1.08E-08 -4.50E-07  9.20E-08  2.78E-06  3.85E-06 -6.79E-07 -7.67E-07  5.50E-07  1.04E-06 -2.10E-05  2.05E-07
          2.14E-07  8.12E-07  3.13E-07  1.81E-07  7.20E-07  2.95E-07  3.98E-06
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        2.81E-06  1.80E-06 -6.14E-07 -5.17E-07  1.38E-06  6.03E-06 -3.20E-07 -1.59E-06  7.08E-07  6.19E-07 -8.85E-06  1.96E-07
          1.71E-07  5.32E-07  4.62E-07  1.16E-07  5.76E-07  3.86E-07  2.05E-06  2.17E-06
 
 OM44
+        2.12E-06  3.41E-06 -3.42E-07 -9.53E-07  1.38E-06  7.50E-06 -3.91E-07 -2.03E-06  8.06E-07  6.83E-07 -1.63E-06  1.72E-07
          1.19E-07  3.02E-07  4.65E-07  1.29E-07  2.90E-07  5.40E-07  9.64E-07  1.74E-06  2.71E-06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 ********************              CORRELATION MATRIX OF ESTIMATE (From Sample Variance)             ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        3.21E-02
 
 TH 2
+        1.45E-02  2.86E-02
 
 TH 3
+       -9.67E-01 -7.24E-03  9.36E-03
 
 TH 4
+       -5.61E-03 -9.58E-01  2.05E-03  8.32E-03
 
 TH 5
+        2.02E-01  2.05E-02 -1.81E-01 -2.75E-03  3.81E-02
 
 TH 6
+        2.40E-02  2.60E-01 -1.65E-02 -2.29E-01  2.52E-02  3.47E-02
 
 TH 7
+       -1.92E-01 -2.24E-02  1.82E-01  4.28E-03 -9.67E-01 -3.02E-02  1.10E-02
 
 TH 8
+       -1.59E-02 -2.44E-01  1.06E-02  2.39E-01 -1.15E-02 -9.58E-01  1.58E-02  1.01E-02
 
 TH 9
+        1.07E-01  1.34E-01 -8.20E-02 -3.08E-02  1.53E-01  2.19E-01 -1.55E-01 -1.46E-01  1.01E-02
 
 TH10
+        8.72E-02  1.16E-01 -6.03E-02 -3.34E-02  1.37E-01  2.21E-01 -1.26E-01 -1.38E-01  6.63E-01  8.33E-03
 
 TH11
+        1.85E-02  2.60E-02 -1.89E-02 -2.00E-02 -4.96E-02 -2.19E-02  4.92E-02  3.60E-02 -1.69E-04 -4.81E-02  9.51E-02
 
 OM11
+        8.59E-04 -2.99E-02  3.28E-03  3.30E-02  2.34E-02 -7.33E-03 -1.39E-02  8.23E-03  4.66E-02 -7.33E-03  9.89E-03  9.47E-04
 
 OM12
+        3.26E-02  1.33E-03 -2.88E-02  3.21E-03  5.51E-02  5.14E-03 -4.50E-02  1.48E-02  3.55E-02  4.91E-02  8.89E-04  3.20E-01
          7.96E-04
 
 OM13
+        1.17E-02 -5.23E-03  2.13E-03 -9.72E-04  1.93E-02 -1.97E-02 -2.02E-03  2.38E-02  3.81E-02  5.12E-04 -7.04E-02  4.05E-01
          2.00E-01  1.10E-03
 
 OM14
+        4.33E-02  2.49E-03 -2.80E-02  2.60E-04  3.24E-02  1.44E-03 -1.58E-02 -2.15E-03  6.26E-02  1.17E-02 -4.43E-02  3.39E-01
          2.21E-01  6.28E-01  9.19E-04
 
 OM22
+        8.96E-04  2.31E-02  4.09E-04 -1.04E-02  2.51E-02 -1.13E-02 -1.95E-02  2.92E-02  3.51E-02  3.17E-02  6.93E-03  1.28E-02
          2.57E-01 -2.82E-02  2.04E-02  1.30E-03
 
 OM23
+        1.30E-02 -5.12E-02 -5.42E-03  4.39E-02  1.25E-02 -2.92E-02  4.30E-03  3.75E-02  5.43E-02  9.00E-03  1.50E-02  1.07E-01
          2.94E-01  3.04E-01  1.65E-01  1.30E-01  1.31E-03
 
 OM24
+        2.03E-02 -2.27E-02 -8.08E-03  1.95E-02  4.91E-02 -1.76E-02 -4.57E-02  2.58E-02  7.90E-02  3.43E-02  2.41E-02  9.08E-02
          1.81E-01  1.72E-01  2.62E-01  2.15E-01  5.76E-01  1.05E-03
 
 OM33
+        3.45E-02  1.90E-04 -2.41E-02  5.55E-03  3.65E-02  5.57E-02 -3.09E-02 -3.79E-02  2.72E-02  6.25E-02 -1.11E-01  1.08E-01
          1.35E-01  3.70E-01  1.70E-01  6.97E-02  2.76E-01  1.41E-01  2.00E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        5.93E-02  4.28E-02 -4.45E-02 -4.21E-02  2.46E-02  1.18E-01 -1.97E-02 -1.06E-01  4.75E-02  5.04E-02 -6.31E-02  1.41E-01
          1.46E-01  3.28E-01  3.41E-01  6.07E-02  2.99E-01  2.49E-01  6.96E-01  1.47E-03
 
 OM44
+        4.00E-02  7.23E-02 -2.22E-02 -6.96E-02  2.19E-02  1.31E-01 -2.16E-02 -1.22E-01  4.84E-02  4.98E-02 -1.04E-02  1.10E-01
          9.06E-02  1.66E-01  3.07E-01  6.02E-02  1.34E-01  3.11E-01  2.93E-01  7.15E-01  1.65E-03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                              MCMC BAYESIAN ANALYSIS                            ********************
 ********************           INVERSE COVARIANCE MATRIX OF ESTIMATE (From Sample Variance)         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.56E+04
 
 TH 2
+       -2.57E+02  1.83E+04
 
 TH 3
+        5.16E+04 -1.06E+03  1.82E+05
 
 TH 4
+       -7.35E+02  6.01E+04 -3.10E+03  2.13E+05
 
 TH 5
+       -2.32E+03  1.39E+01 -7.43E+03  3.50E+01  1.12E+04
 
 TH 6
+        7.41E+01 -3.51E+03  2.18E+02 -1.15E+04  2.43E+02  1.21E+04
 
 TH 7
+       -7.23E+03  1.31E+02 -2.48E+04  3.72E+02  3.74E+04  8.45E+02  1.34E+05
 
 TH 8
+        2.46E+02 -1.14E+04  7.64E+02 -4.03E+04  6.25E+02  3.90E+04  2.29E+03  1.36E+05
 
 TH 9
+       -7.02E+02 -4.15E+03 -1.67E+03 -1.37E+04  7.25E+02 -7.02E+02  3.53E+03 -1.15E+03  1.93E+04
 
 TH10
+       -8.41E+02 -6.11E+02 -2.72E+03 -1.95E+03 -1.25E+03 -2.96E+03 -3.87E+03 -8.74E+03 -1.38E+04  2.74E+04
 
 TH11
+       -1.55E+01 -2.82E+01 -1.48E+01 -5.15E+01  2.56E+00 -6.45E+01 -4.37E+01 -2.49E+02 -4.57E+01  1.15E+02  1.14E+02
 
 OM11
+        2.20E+03 -9.09E+02  5.03E+03 -9.77E+03 -1.70E+03  2.54E+03 -5.37E+03  9.17E+03 -9.21E+03  9.38E+03 -4.85E+02  1.48E+06
 
 OM12
+        4.40E+02  1.24E+03  4.80E+03  7.18E+03 -7.17E+02 -7.58E+03  2.80E+03 -2.52E+04  9.69E+03 -1.24E+04 -3.39E+00 -4.83E+05
          2.08E+06
 
 OM13
+       -2.77E+03  6.30E+03 -1.27E+04  2.79E+04 -1.41E+03 -2.63E+03 -5.58E+03 -1.31E+04 -1.89E+03  3.10E+03  5.14E+02 -4.29E+05
          1.08E+05  1.82E+06
 
 OM14
+       -2.59E+03 -5.94E+03 -4.26E+03 -2.40E+04 -6.21E+03  5.57E+03 -2.20E+04  1.66E+04 -7.12E+03  5.40E+03  3.45E+02 -1.11E+05
         -2.64E+05 -1.20E+06  2.39E+06
 
 OM22
+        1.02E+03 -3.63E+03  2.35E+03 -9.76E+03 -1.96E+03 -2.11E+03 -6.88E+03 -9.54E+03 -1.17E+03  1.45E+03  1.42E+01  4.31E+04
         -2.86E+05  7.67E+04  1.39E+04  6.71E+05
 
 OM23
+        1.25E+03  3.17E+03  3.60E+03  5.34E+03 -6.44E+03  3.30E+02 -2.60E+04 -3.41E+02 -6.96E+03  6.80E+03 -3.67E+02  9.65E+04
         -3.51E+05 -3.47E+05  3.16E+05  3.26E+04  1.12E+06
 
 OM24
+       -2.81E+03  5.39E+03 -1.06E+04  1.70E+04  4.81E+03  6.81E+02  2.21E+04 -2.27E+03 -5.96E+03 -5.37E+02 -2.04E+02 -2.21E+04
          1.23E+05  2.13E+05 -4.20E+05 -1.78E+05 -7.76E+05  1.65E+06
 
 OM33
+       -8.92E+02 -3.23E+03 -3.84E+03 -1.36E+04  4.62E+01 -2.33E+03  2.19E+03 -8.62E+03  5.60E+03 -7.30E+03  6.07E+02  3.39E+04
         -3.41E+04 -3.01E+05  2.80E+05 -3.45E+04  2.75E+04 -4.72E+04  6.50E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        4.08E+02  4.83E+03  7.58E+03  1.85E+04  1.34E+03 -2.09E+02  2.39E+03  5.83E+03 -4.13E+03  5.97E+03 -1.07E+01 -8.60E+03
          4.22E+04  1.56E+05 -4.26E+05  1.50E+04 -3.35E+05  2.76E+05 -8.24E+05  2.18E+06
 
 OM44
+       -4.08E+03 -2.13E+03 -1.59E+04 -3.59E+03  7.75E+02 -2.59E+03  3.02E+03 -3.72E+03  3.89E+03 -3.87E+03 -1.51E+02 -2.35E+04
         -1.11E+04  2.50E+04 -3.67E+04  9.70E+02  2.28E+05 -3.53E+05  2.90E+05 -1.06E+06  1.00E+06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 
 
 #TBLN:      6
 #METH: First Order Conditional Estimation with Interaction (No Prior)
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
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
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      14
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     14
 NOPRIOR SETTING (NOPRIOR):                 ON
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          ON
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      OFF
 RAW OUTPUT FILE (FILE): example2.ext
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

 
0ITERATION NO.:    0    OBJECTIVE VALUE:  -7981.93164712434        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  3.3009E+00  3.2532E+00 -6.1243E-01 -2.0833E-01  7.3503E-01  1.1396E+00  3.3504E-01  1.9158E-01  6.9166E-01  2.2992E+00
            -3.7204E-02  1.0296E-02  3.7205E-06  9.2078E-04 -8.9716E-04  7.7895E-03 -5.2439E-04  2.7619E-04  8.4949E-03  7.7703E-04
             8.9914E-03
 PARAMETER:  1.0000E-01  1.0000E-01 -1.0000E-01 -1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
            -1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01 -1.0000E-01  1.0000E-01 -1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.9322E+03 -2.1558E+03 -1.5704E+03 -3.3011E+02  1.6673E+02 -8.2229E+00 -3.2705E+02  2.8527E+02  2.0318E+03  2.7926E+03
             7.1750E+04 -1.1839E+02 -6.7168E-01 -8.2885E+01 -1.0965E+02 -6.4461E+02 -1.6193E+01 -1.2515E+01 -7.7275E+02 -4.3678E+02
            -4.1072E+02
 
0ITERATION NO.:    5    OBJECTIVE VALUE:  -10691.3255247918        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:       97
 NPARAMETR:  3.3289E+00  3.2831E+00 -6.0495E-01 -2.0836E-01  7.3087E-01  1.1299E+00  3.3693E-01  1.8979E-01  6.8428E-01  2.2756E+00
            -9.3321E-02  1.0295E-02  3.7204E-06  9.2198E-04 -8.9373E-04  7.8108E-03 -5.2494E-04  2.7663E-04  8.5224E-03  7.8681E-04
             9.0065E-03
 PARAMETER:  1.0085E-01  1.0092E-01 -9.8779E-02 -1.0002E-01  9.9433E-02  9.9149E-02  1.0057E-01  9.9062E-02  9.8934E-02  9.8975E-02
            -2.5084E-01  9.9954E-02  1.0000E-01  1.0014E-01 -9.9622E-02  1.0137E-01 -9.9968E-02  1.0002E-01  1.0163E-01  1.0093E-01
             1.0079E-01
 GRADIENT:   5.2393E+04  3.7048E+04  3.3047E+04  7.7749E+03 -1.4396E+03 -4.2746E+03 -2.3679E+03 -2.4813E+03  4.0690E+02 -2.6153E+04
             1.2235E+03 -1.0507E+02  3.7224E-02  3.6509E-01  3.1718E+01 -3.3435E+01  4.4305E+00 -2.0404E+00 -1.4502E+01 -1.3772E+01
            -1.9777E+01
 
0ITERATION NO.:   10    OBJECTIVE VALUE:  -10746.2517333300        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      177
 NPARAMETR:  3.2919E+00  3.2607E+00 -6.0400E-01 -2.0693E-01  7.6344E-01  1.1259E+00  3.2542E-01  1.9631E-01  6.9017E-01  2.2941E+00
            -1.0432E-01  7.7523E-03  3.2297E-06  7.1342E-04 -2.0662E-04  7.8254E-03 -5.2420E-04  2.7226E-04  8.5109E-03  8.8814E-04
             8.2757E-03
 PARAMETER:  9.9729E-02  1.0023E-01 -9.8623E-02 -9.9332E-02  1.0387E-01  9.8800E-02  9.7128E-02  1.0247E-01  9.9785E-02  9.9779E-02
            -2.8039E-01 -4.1865E-02  1.0004E-01  8.9290E-02 -2.6540E-02  1.0230E-01 -9.9729E-02  9.8265E-02  1.0196E-01  1.0545E-01
             6.1506E-02
 GRADIENT:   1.9020E+04  1.1994E+04  1.3382E+04  2.4604E+03 -1.7428E+03  1.9168E+02 -3.4036E+03  2.3211E+02 -3.7537E+02 -8.9838E+03
            -3.1289E+02 -1.8708E+02  1.1052E-01 -6.9370E+00  6.5998E+01  5.2086E+00  1.4287E+00 -4.8621E+00 -1.1780E-01  2.3566E+00
            -2.9250E+01
 
0ITERATION NO.:   15    OBJECTIVE VALUE:  -10747.6725605324        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  3.2874E+00  3.2577E+00 -6.0302E-01 -2.0598E-01  7.3244E-01  1.1470E+00  3.3438E-01  1.9007E-01  6.9151E-01  2.2969E+00
            -1.0391E-01  8.2439E-03  3.3263E-06  1.0194E-03  5.3322E-04  7.8915E-03 -5.6147E-04  3.0396E-04  9.0560E-03  7.6485E-04
             8.4920E-03
 PARAMETER:  9.9593E-02  1.0014E-01 -9.8464E-02 -9.8875E-02  9.9647E-02  1.0065E-01  9.9804E-02  9.9208E-02  9.9980E-02  9.9898E-02
            -2.7930E-01 -1.1124E-02  9.9913E-02  1.2373E-01  6.6420E-02  1.0651E-01 -1.0639E-01  1.0914E-01  1.2970E-01  7.9856E-02
             7.5318E-02
 GRADIENT:   1.6929E+04  1.1075E+04  1.2029E+04  2.3781E+03 -1.8837E+03  2.2285E+02 -3.1205E+03  1.4198E+01 -5.0590E+02 -6.6951E+03
            -3.2349E+02 -1.6964E+02  6.6836E-02 -1.8277E+01  9.1963E+01  8.3079E+00  1.8205E+00 -5.0841E+00  1.6671E+01 -1.4063E+01
            -2.4228E+01
 
0ITERATION NO.:   20    OBJECTIVE VALUE:  -10754.5396412447        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  3.3009E+00  3.2591E+00 -6.0727E-01 -2.0643E-01  7.3584E-01  1.1466E+00  3.3463E-01  1.9038E-01  6.9630E-01  2.2994E+00
            -1.0126E-01  8.9930E-03  3.5996E-06  2.9063E-03  1.4851E-03  7.0646E-03  1.1951E-03  1.2929E-03  1.5214E-02  5.8716E-03
             1.2583E-02
 PARAMETER:  1.0000E-01  1.0018E-01 -9.9158E-02 -9.9090E-02  1.0011E-01  1.0062E-01  9.9876E-02  9.9372E-02  1.0067E-01  1.0001E-01
            -2.7219E-01  3.2365E-02  1.0352E-01  3.3772E-01  1.7712E-01  5.1159E-02  2.3893E-01  4.9074E-01  3.5938E-01  4.5569E-01
             1.7311E-01
 GRADIENT:   1.3774E+04  8.8493E+03  9.3101E+03  1.9200E+03 -1.1382E+03  5.5446E+01 -1.8596E+03 -4.8518E+00  1.6769E+02 -6.3003E+03
            -1.9560E+02 -1.2717E+02  1.0149E-02 -2.7335E+00  6.2231E+01 -1.0675E+01  5.5126E+00 -2.1050E+00  2.3589E+01 -6.8380E+00
            -1.8719E+00
 
0ITERATION NO.:   25    OBJECTIVE VALUE:  -10772.0656579881        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      395
 NPARAMETR:  3.3053E+00  3.2580E+00 -6.1186E-01 -2.0825E-01  7.3515E-01  1.1398E+00  3.3621E-01  1.9177E-01  6.9429E-01  2.3013E+00
            -1.0035E-01  1.0191E-02  7.9038E-06  1.1884E-03 -6.5352E-04  7.8515E-03 -3.2682E-04  4.6340E-04  9.5618E-03  1.6737E-03
             9.3807E-03
 PARAMETER:  1.0013E-01  1.0015E-01 -9.9908E-02 -9.9965E-02  1.0002E-01  1.0002E-01  1.0035E-01  1.0010E-01  1.0038E-01  1.0009E-01
            -2.6974E-01  9.4916E-02  2.1352E-01  1.2972E-01 -7.3215E-02  1.0396E-01 -6.2213E-02  1.6711E-01  1.5811E-01  1.9060E-01
             1.0938E-01
 GRADIENT:   7.4996E+01  6.5226E+01  4.2213E+01  1.1240E+01 -4.5002E+01  1.6262E+01 -7.1100E+01  5.1751E+00 -9.1320E+00 -1.0931E+02
            -6.3224E+00 -1.2441E+00 -1.9100E-02 -3.9459E-02  1.4903E+00  3.1306E-01 -1.5774E-01  8.4809E-02 -1.2442E+00  6.7218E-02
            -2.8055E+00
 
0ITERATION NO.:   30    OBJECTIVE VALUE:  -10772.1218018054        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  3.3058E+00  3.2583E+00 -6.1192E-01 -2.0826E-01  7.3518E-01  1.1403E+00  3.3627E-01  1.9170E-01  6.9493E-01  2.3017E+00
            -1.0016E-01  1.0252E-02  1.0227E-04  1.2232E-03 -6.1744E-04  7.9281E-03 -2.1035E-04  5.0820E-04  9.9379E-03  1.9516E-03
             9.5913E-03
 PARAMETER:  1.0015E-01  1.0016E-01 -9.9917E-02 -9.9970E-02  1.0002E-01  1.0006E-01  1.0037E-01  1.0006E-01  1.0047E-01  1.0011E-01
            -2.6923E-01  9.7853E-02  2.7546E+00  1.3313E-01 -6.8970E-02  1.0875E-01 -4.2045E-02  1.8439E-01  1.7771E-01  2.1546E-01
             1.1585E-01
 GRADIENT:   2.8101E+02  1.1342E+02  1.6449E+02  2.1795E+01 -3.4396E+01 -6.5305E+00 -4.8765E+01 -5.2130E+00  5.6710E+01 -1.0415E+02
            -9.2929E+00  6.4131E-01 -1.0771E-02 -1.2385E+00  1.2837E+00  7.6649E-01  1.8158E-01 -1.2068E-01  4.0597E-02  3.7089E-01
            -1.6112E+00
 
0ITERATION NO.:   35    OBJECTIVE VALUE:  -10772.1430212060        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  3.3052E+00  3.2580E+00 -6.1182E-01 -2.0823E-01  7.3532E-01  1.1399E+00  3.3630E-01  1.9180E-01  6.9454E-01  2.3016E+00
            -1.0003E-01  1.0282E-02  1.9104E-04  1.2814E-03 -5.9571E-04  7.9490E-03 -1.7050E-04  5.4236E-04  1.0007E-02  2.0014E-03
             9.6711E-03
 PARAMETER:  1.0013E-01  1.0015E-01 -9.9900E-02 -9.9952E-02  1.0004E-01  1.0003E-01  1.0037E-01  1.0011E-01  1.0042E-01  1.0010E-01
            -2.6887E-01  9.9324E-02  5.1381E+00  1.3926E-01 -6.6444E-02  1.0991E-01 -3.6666E-02  1.9817E-01  1.8058E-01  2.2005E-01
             1.1907E-01
 GRADIENT:  -7.9383E+00 -9.6817E+00 -3.7708E+00 -1.1959E+00 -1.3189E-01 -2.1198E-02 -2.8028E-01  3.7267E-03 -3.1691E-01 -1.0338E+01
            -2.1318E-02  5.8793E-03 -1.1027E-04 -8.0300E-03  3.8584E-03  3.3304E-03  2.9661E-03 -1.1804E-03  9.0481E-03  1.2014E-03
            -1.7269E-02
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      566
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         7.0199E-05 -1.0238E-03 -3.6824E-05 -9.6585E-04
 SE:             4.7162E-03  2.9847E-03  2.9483E-03  3.6726E-03
 N:                     400         400         400         400
 
 P VAL.:         9.8812E-01  7.3158E-01  9.9003E-01  7.9256E-01
 
 ETASHRINKSD(%)  6.8614E+00  3.2963E+01  4.0980E+01  2.5216E+01
 ETASHRINKVR(%)  1.3252E+01  5.5061E+01  6.5166E+01  4.4074E+01
 EBVSHRINKSD(%)  6.9635E+00  3.3079E+01  4.1033E+01  2.5321E+01
 EBVSHRINKVR(%)  1.3442E+01  5.5216E+01  6.5229E+01  4.4230E+01
 EPSSHRINKSD(%)  2.6153E+01
 EPSSHRINKVR(%)  4.5466E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2000
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    3675.75413281869     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -10772.1430212060     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -7096.38888838734     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                          1600
  
 #TERE:
 Elapsed estimation  time in seconds:    18.40
 Elapsed covariance  time in seconds:    21.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************          FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR)        ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   -10772.143       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************          FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR)        ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.31E+00  3.26E+00 -6.12E-01 -2.08E-01  7.35E-01  1.14E+00  3.36E-01  1.92E-01  6.95E-01  2.30E+00 -1.00E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.03E-02
 
 ETA2
+        1.91E-04  7.95E-03
 
 ETA3
+        1.28E-03 -1.71E-04  1.00E-02
 
 ETA4
+       -5.96E-04  5.42E-04  2.00E-03  9.67E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        1.01E-01
 
 ETA2
+        2.11E-02  8.92E-02
 
 ETA3
+        1.26E-01 -1.91E-02  1.00E-01
 
 ETA4
+       -5.97E-02  6.19E-02  2.03E-01  9.83E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************          FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR)        ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.27E-02  2.86E-02  9.52E-03  8.30E-03  3.92E-02  3.58E-02  1.13E-02  1.04E-02  1.05E-02  8.57E-03  2.83E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        9.68E-04
 
 ETA2
+        8.28E-04  1.38E-03
 
 ETA3
+        1.28E-03  1.44E-03  3.05E-03
 
 ETA4
+        1.00E-03  1.10E-03  2.15E-03  1.93E-03
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4     
 
 ETA1
+        4.77E-03
 
 ETA2
+        9.07E-02  7.71E-03
 
 ETA3
+        1.15E-01  1.62E-01  1.52E-02
 
 ETA4
+        1.04E-01  1.22E-01  1.76E-01  9.80E-03
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************          FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR)        ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.07E-03
 
 TH 2
+        1.55E-05  8.19E-04
 
 TH 3
+       -3.01E-04 -2.84E-06  9.06E-05
 
 TH 4
+       -9.82E-07 -2.28E-04  1.00E-07  6.89E-05
 
 TH 5
+        2.90E-04  2.48E-05 -7.84E-05 -2.42E-06  1.53E-03
 
 TH 6
+        2.82E-05  2.92E-04 -5.49E-06 -7.60E-05  3.98E-05  1.28E-03
 
 TH 7
+       -8.09E-05 -7.82E-06  2.28E-05  8.15E-07 -4.30E-04 -1.24E-05  1.28E-04
 
 TH 8
+       -4.56E-06 -7.80E-05  9.16E-07  2.23E-05 -5.45E-06 -3.55E-04  1.85E-06  1.07E-04
 
 TH 9
+        4.38E-05  4.51E-05 -9.15E-06 -4.54E-06  5.74E-05  7.05E-05 -1.80E-05 -1.18E-05  1.10E-04
 
 TH10
+        2.72E-05  3.16E-05 -4.97E-06 -3.50E-06  5.17E-05  5.99E-05 -1.50E-05 -9.26E-06  6.31E-05  7.35E-05
 
 TH11
+       -1.38E-06  2.58E-07  3.27E-07 -1.04E-07 -3.30E-06  2.83E-06  1.07E-06 -7.06E-07 -5.06E-07 -8.65E-08  7.99E-06
 
 OM11
+        2.65E-07 -1.01E-07 -6.72E-09  7.12E-08  3.48E-07  4.68E-07  2.14E-09 -3.01E-08  4.27E-07  3.99E-07  4.00E-07  9.37E-07
 
 OM12
+        2.18E-07 -1.97E-07  3.81E-09  1.23E-07  1.36E-07 -1.92E-07  1.46E-07  2.03E-07  3.41E-07  3.21E-07  4.34E-07  2.74E-07
          6.86E-07
 
 OM13
+        4.80E-07  1.34E-07  4.31E-08  3.73E-08  1.12E-06  1.17E-06 -1.61E-07 -1.22E-07  1.01E-06  1.00E-06  7.31E-07  5.41E-07
          2.18E-07  1.63E-06
 
 OM14
+        9.59E-07 -1.32E-06 -1.93E-07  5.03E-07  6.24E-07  9.61E-07 -8.91E-08 -1.19E-07  8.93E-07  6.87E-07  6.13E-07  3.63E-07
          2.11E-07  8.99E-07  1.01E-06
 
 OM22
+       -2.57E-07 -9.91E-08  1.50E-07  1.13E-07 -1.28E-06 -1.20E-06  7.16E-07  5.50E-07  1.29E-07  1.68E-07  1.24E-06  8.56E-08
          3.70E-07  7.84E-08  7.94E-08  1.89E-06
 
 OM23
+        6.60E-07  1.01E-07 -7.04E-08  8.09E-08  2.05E-06 -5.83E-07 -2.40E-07  3.87E-07  1.20E-06  9.22E-07  5.82E-07  1.75E-07
          3.54E-07  6.75E-07  3.60E-07  3.83E-07  2.06E-06
 
 OM24
+        2.88E-07 -2.84E-07 -1.34E-08  1.45E-07  1.67E-06 -2.49E-06 -3.18E-07  8.97E-07  6.30E-07  4.78E-07  3.95E-07  9.44E-08
          2.31E-07  2.88E-07  3.22E-07  4.17E-07  9.98E-07  1.20E-06
 
 OM33
+        7.59E-07  1.73E-06  5.55E-08 -1.12E-07  1.34E-06  5.42E-06 -6.84E-08 -9.47E-07  1.53E-06  1.93E-06  3.50E-06  4.56E-07
          3.37E-07  1.96E-06  1.08E-06  3.86E-07  1.62E-06  6.31E-07  9.29E-06
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        4.55E-07  2.11E-06  6.03E-08 -3.10E-07  3.77E-07  5.77E-06  4.97E-08 -1.23E-06  1.68E-06  1.53E-06  2.50E-06  3.05E-07
          2.42E-07  1.23E-06  9.37E-07  2.27E-07  1.16E-06  6.22E-07  5.58E-06  4.60E-06
 
 OM44
+        5.39E-08  2.05E-06  1.29E-07 -3.71E-07 -3.94E-07  5.71E-06  1.74E-07 -1.35E-06  1.33E-06  1.20E-06  1.91E-06  2.03E-07
          1.68E-07  7.23E-07  7.63E-07  1.48E-07  6.69E-07  6.48E-07  3.15E-06  3.33E-06  3.71E-06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************          FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR)        ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        3.27E-02
 
 TH 2
+        1.66E-02  2.86E-02
 
 TH 3
+       -9.68E-01 -1.04E-02  9.52E-03
 
 TH 4
+       -3.62E-03 -9.59E-01  1.27E-03  8.30E-03
 
 TH 5
+        2.27E-01  2.21E-02 -2.10E-01 -7.45E-03  3.92E-02
 
 TH 6
+        2.42E-02  2.85E-01 -1.61E-02 -2.56E-01  2.84E-02  3.58E-02
 
 TH 7
+       -2.19E-01 -2.41E-02  2.11E-01  8.67E-03 -9.69E-01 -3.06E-02  1.13E-02
 
 TH 8
+       -1.35E-02 -2.63E-01  9.29E-03  2.60E-01 -1.34E-02 -9.59E-01  1.58E-02  1.04E-02
 
 TH 9
+        1.28E-01  1.50E-01 -9.18E-02 -5.22E-02  1.40E-01  1.88E-01 -1.52E-01 -1.09E-01  1.05E-02
 
 TH10
+        9.72E-02  1.29E-01 -6.09E-02 -4.91E-02  1.54E-01  1.95E-01 -1.54E-01 -1.04E-01  7.03E-01  8.57E-03
 
 TH11
+       -1.50E-02  3.19E-03  1.22E-02 -4.42E-03 -2.98E-02  2.79E-02  3.34E-02 -2.41E-02 -1.71E-02 -3.57E-03  2.83E-03
 
 OM11
+        8.38E-03 -3.66E-03 -7.29E-04  8.85E-03  9.17E-03  1.35E-02  1.95E-04 -3.00E-03  4.21E-02  4.81E-02  1.46E-01  9.68E-04
 
 OM12
+        8.05E-03 -8.30E-03  4.83E-04  1.79E-02  4.20E-03 -6.48E-03  1.56E-02  2.37E-02  3.93E-02  4.51E-02  1.85E-01  3.41E-01
          8.28E-04
 
 OM13
+        1.15E-02  3.68E-03  3.55E-03  3.52E-03  2.24E-02  2.57E-02 -1.11E-02 -9.23E-03  7.52E-02  9.14E-02  2.03E-01  4.38E-01
          2.07E-01  1.28E-03
 
 OM14
+        2.93E-02 -4.58E-02 -2.02E-02  6.03E-02  1.59E-02  2.67E-02 -7.84E-03 -1.15E-02  8.50E-02  7.98E-02  2.16E-01  3.74E-01
          2.54E-01  7.01E-01  1.00E-03
 
 OM22
+       -5.72E-03 -2.52E-03  1.15E-02  9.85E-03 -2.39E-02 -2.44E-02  4.60E-02  3.86E-02  8.98E-03  1.43E-02  3.20E-01  6.43E-02
          3.25E-01  4.46E-02  5.75E-02  1.38E-03
 
 OM23
+        1.41E-02  2.45E-03 -5.16E-03  6.79E-03  3.65E-02 -1.14E-02 -1.48E-02  2.60E-02  7.98E-02  7.49E-02  1.43E-01  1.26E-01
          2.98E-01  3.68E-01  2.50E-01  1.94E-01  1.44E-03
 
 OM24
+        8.04E-03 -9.03E-03 -1.28E-03  1.60E-02  3.89E-02 -6.34E-02 -2.56E-02  7.89E-02  5.48E-02  5.08E-02  1.27E-01  8.89E-02
          2.54E-01  2.05E-01  2.92E-01  2.76E-01  6.34E-01  1.10E-03
 
 OM33
+        7.63E-03  1.98E-02  1.91E-03 -4.41E-03  1.12E-02  4.97E-02 -1.98E-03 -3.00E-02  4.80E-02  7.39E-02  4.06E-01  1.55E-01
          1.34E-01  5.04E-01  3.53E-01  9.21E-02  3.71E-01  1.89E-01  3.05E-03
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        6.49E-03  3.44E-02  2.95E-03 -1.74E-02  4.49E-03  7.52E-02  2.05E-03 -5.54E-02  7.50E-02  8.30E-02  4.12E-01  1.47E-01
          1.36E-01  4.49E-01  4.35E-01  7.68E-02  3.78E-01  2.64E-01  8.53E-01  2.15E-03
 
 OM44
+        8.56E-04  3.71E-02  7.06E-03 -2.32E-02 -5.22E-03  8.27E-02  7.96E-03 -6.75E-02  6.57E-02  7.26E-02  3.50E-01  1.09E-01
          1.06E-01  2.94E-01  3.94E-01  5.59E-02  2.42E-01  3.07E-01  5.36E-01  8.05E-01  1.93E-03
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************          FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (NO PRIOR)        ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 TH 1
+        1.59E+04
 
 TH 2
+        0.00E+00  1.85E+04
 
 TH 3
+        5.24E+04  0.00E+00  1.84E+05
 
 TH 4
+        0.00E+00  6.08E+04  0.00E+00  2.15E+05
 
 TH 5
+       -2.14E+03  0.00E+00 -6.77E+03  0.00E+00  1.12E+04
 
 TH 6
+        2.41E-06 -3.51E+03  1.69E-21 -1.13E+04  1.06E-22  1.17E+04
 
 TH 7
+       -6.77E+03  0.00E+00 -2.34E+04  0.00E+00  3.76E+04  0.00E+00  1.34E+05
 
 TH 8
+        0.00E+00 -1.13E+04  0.00E+00 -3.94E+04  0.00E+00  3.83E+04  0.00E+00  1.35E+05
 
 TH 9
+       -1.36E+03 -3.88E+03 -3.62E+03 -1.25E+04  1.09E+03 -5.66E+02  4.23E+03 -1.28E+03  1.95E+04
 
 TH10
+       -9.76E+02 -2.73E+02 -3.56E+03 -8.22E+02 -7.37E+02 -3.20E+03 -1.43E+03 -9.86E+03 -1.49E+04  2.83E+04
 
 TH11
+        1.46E+03  1.36E+03  4.72E+03  4.74E+03  1.11E+03  4.95E+02  3.40E+03  2.26E+03  1.09E+03  6.94E+02  1.76E+05
 
 OM11
+       -7.04E+01 -2.46E+02 -1.65E+02 -7.69E+02 -1.36E+03 -5.10E+02 -4.67E+03 -1.41E+03 -8.13E+02 -4.27E+02 -3.48E+04  1.48E+06
 
 OM12
+       -1.89E+03 -1.03E+03 -6.65E+03 -3.07E+03 -4.34E+03 -3.68E+03 -1.57E+04 -1.36E+04  1.02E+03 -1.63E+03 -2.18E+04 -4.87E+05
          1.99E+06
 
 OM13
+       -5.14E+03  5.45E+03 -2.12E+04  2.77E+04  2.88E+02 -9.57E+02  1.79E+03 -6.02E+03  1.58E+03 -5.60E+03  2.71E+04 -4.50E+05
          1.71E+05  1.73E+06
 
 OM14
+        2.63E+03 -3.47E+03  1.47E+04 -3.07E+04 -1.19E+03  7.40E+02 -5.37E+03  6.85E+03 -6.22E+03  3.26E+03 -3.14E+04 -1.00E+05
         -2.97E+05 -1.28E+06  2.44E+06
 
 OM22
+       -1.21E+03 -2.00E+03 -4.12E+03 -6.90E+03 -5.32E+03 -2.22E+03 -2.13E+04 -8.45E+03 -4.04E+02 -5.65E+01 -1.05E+05  4.64E+04
         -2.95E+05  8.63E+02  7.24E+04  6.94E+05
 
 OM23
+        1.61E+03  3.49E+02  6.42E+03 -1.51E+03 -5.92E+03  1.42E+03 -2.09E+04  7.61E+03 -4.88E+03  3.34E+02  1.40E+04  8.01E+04
         -2.62E+05 -3.57E+05  3.23E+05  2.28E+04  1.06E+06
 
 OM24
+       -6.04E+02  1.65E+03 -3.82E+03  9.70E+03  1.87E+03 -7.08E+02  1.08E+04 -1.47E+04  6.66E+02 -2.92E+02  2.51E+04  1.19E+04
          7.34E+03  2.62E+05 -4.59E+05 -2.18E+05 -8.16E+05  1.70E+06
 
 OM33
+       -5.95E+02 -2.29E+03 -1.30E+03 -1.02E+04 -1.02E+03 -1.04E+03 -3.04E+03 -3.93E+03  5.23E+03 -4.08E+03 -4.80E+04  3.98E+04
         -1.64E+04 -2.97E+05  2.36E+05 -1.38E+03  1.57E+04 -1.40E+04  5.87E+05
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM22      OM23      OM24      OM33      OM34      OM44      SG11  
 
 OM34
+        4.02E+02  3.13E+02  5.55E+02  5.33E+03  1.53E+03 -7.25E+02  5.04E+03 -2.19E+03 -8.06E+03  6.39E+03 -5.89E+03  1.18E+04
          5.93E+04  2.13E+05 -3.85E+05  1.49E+04 -2.97E+05  2.42E+05 -9.11E+05  2.25E+06
 
 OM44
+       -1.36E+03 -9.54E+02 -5.10E+03 -2.86E+03 -1.16E+02 -6.17E+02 -1.68E+03  2.31E+03  3.20E+03 -4.39E+03 -4.40E+04  5.56E+03
         -4.70E+03  9.95E+03 -5.09E+04  4.55E+04  2.08E+05 -3.20E+05  3.51E+05 -1.19E+06  1.09E+06
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.48
 #CPUT: Total CPU Time in Seconds,      558.328
Stop Time: 
Sat 04/22/2017 
09:26 AM
