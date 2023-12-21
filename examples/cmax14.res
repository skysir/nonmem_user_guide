Fri 05/18/2018 
03:35 PM


$PROB  THEOPHYLLINE   POPULATION DATA
$INPUT      ID DOSE=AMT TIME CP=DV WT
$DATA       THEOPPCMAX
$ABB COMRES=2
$ABBR DECLARE INTEGER I

$SUBROUTINES  ADVAN14  TOL=6 OTHER=cmax14u.f90
$MODEL COMP=(DEPOT,INITIALOFF,DEFDOSE) COMP=(CENTRAL,DEFOBS,NOOFF)

$PK

;THETA(1)=MEAN ABSORPTION RATE CONSTANT (1/HR)
;THETA(2)=MEAN ELIMINATION RATE CONSTANT (1/HR)
;THETA(3)=SLOPE OF CLEARANCE VS WEIGHT RELATIONSHIP (LITERS/HR/KG)
;SCALING PARAMETER=VOLUME/WT SINCE DOSE IS WEIGHT-ADJUSTED
   IF (Wt>0) WEIGHT=WT
   KA=THETA(1)+ETA(1)
   KE=THETA(2)+ETA(2)
   CL=THETA(3)*WEIGHT+ETA(3)
   S2=CL/KE/WEIGHT
 

$THETA  (0.0,1.96) (.0,.086) (0.0,.040)
$OMEGA BLOCK(3)  1.519 .00376 .000133 -0.066 .00752 .481

$DES

   DADT(1)=-KA*A(1)
   DADT(2)= KA*A(1)-KE*A(2)


$ERROR
include nonmem_reserved_general

TMAX=0.0
CMAX=0.0
ZRCROSSDIR=0 ; ZERO CROSSING DIRECTION
CMAXDER=0.0
; THIS ROOT-FINDER CAN FIND THE EXACT TMAX AND CMAX, BY FINDING WHEN DADT(2) IS ZERO
; SEE SUBROUTINE FCVROOTFN in cmax14u.f90
 IF(RTSIGNAL>0.0) THEN ; A ZERO-CROSSING HAS BEEN FOUND
  TMAX=TRTINFO(1) ; ONLY ONE CROSSING EXPECTED, SO INDEX IS 1
  CMAX=YRTINFO(2,1)/S2 ; PICK UP THE Y(2) VALUE AT CROSSING NUMBER 1
  CMAXDER=YPRTINFO(2,1); VERIFY THAT CMAX DERIVATIVE IS CLOSE TO ZERO.
  ZRCROSSDIR=RTINFO(1,1) ; VERIFY THAT ZEROCROSSING IS IN NEGATIVE DIRECTION (-1, FROM POSITIVE TO NEGATIVE), INDICATING A PEAK, RATHER THAN A NADIR
 RTSIGNAL=0
 ENDIF

; Compare with analytical solution for tmax, cmax
CMAXC=0.0
TMAXC=0.0
IF(AMT>0.0) THEN
  TMAXC=LOG(KA/KE)/(KA-KE)
  CMAXC=KA*AMT/S2*(EXP(-KE*TMAXC)-EXP(-KA*TMAXC))/(KA-KE)
ENDIF


IPRED=F
   Y=F+EPS(1)

$SIGMA  .476

$EST METHOD=1 MAXEVAL=9999 PRINT=1 SIGL=6 NOABORT
     NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST
$COV MATRIX=R UNCONDITIONAL
$TABLE          ID TIME WEIGHT TMAX CMAX CMAXDER ZRCROSSDIR TMAXC CMAXC IPRED NOAPPEND NOPRINT FORMAT=S1PE15.8 FILE=cmax14.tab
  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
  
License Registered to: IDS NONMEM 7 TEAM
Expiration Date:     2 JUN 2030
Current Date:       18 MAY 2018
Days until program expires :4394
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0 alpha version 2
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 THEOPHYLLINE   POPULATION DATA
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      120
 NO. OF DATA ITEMS IN DATA SET:   7
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   4
 MDV DATA ITEM IS DATA ITEM NO.:  7
0INDICES PASSED TO SUBROUTINE PRED:
   6   3   2   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID DOSE TIME CP WT EVID MDV
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 WEIGHT TMAX CMAX ZRCROSSDIR CMAXDER CMAXC TMAXC IPRED
0FORMAT FOR DATA:
 (5E6.0,2F2.0)

 TOT. NO. OF OBS RECS:      108
 TOT. NO. OF INDIVIDUALS:       12
0LENGTH OF THETA:   3
0DEFAULT THETA BOUNDARY TEST OMITTED:   YES
0OMEGA HAS BLOCK FORM:
  1
  1  1
  1  1  1
0DEFAULT OMEGA BOUNDARY TEST OMITTED:   YES
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:   YES
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1960E+01     0.1000E+07
  0.0000E+00     0.8600E-01     0.1000E+07
  0.0000E+00     0.4000E-01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 BLOCK SET NO.   BLOCK                                                                    FIXED
        1                                                                                   NO
                  0.1519E+01
                  0.3760E-02   0.1330E-03
                 -0.6600E-01   0.7520E-02   0.4810E+00
0INITIAL ESTIMATE OF SIGMA:
 0.4760E+00
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
04 COLUMNS APPENDED:    NO
 PRINTED:                NO
 HEADER:                YES
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE15.8
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID TIME WEIGHT TMAX CMAX CMAXDER ZRCROSSDIR TMAXC CMAXC IPRED
1DOUBLE PRECISION PREDPP VERSION 7.5.0 alpha version 2

 GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF EQUATIONS (CVODES, ADVAN14)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   2
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
 INITIAL (BASE) TOLERANCE SETTINGS:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            3           *           *           *           *
    3            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          3
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   2

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0DES SUBROUTINE USES COMPACT STORAGE MODE.
1
 
 
 #TBLN:      1
 #METH: First Order Conditional Estimation
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     NO
 LAPLACIAN OBJ. FUNC.:                    NO
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
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
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
 RAW OUTPUT FILE (FILE): cmax14.ext
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

 TOLERANCES FOR ESTIMATION/EVALUATION STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR COVARIANCE STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 TOLERANCES FOR TABLE/SCATTER STEP:
 NRD (RELATIVE) VALUE(S) OF TOLERANCE:   6
 ANRD (ABSOLUTE) VALUE(S) OF TOLERANCE:  12
 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=NPRED
 RES=NRES
 WRES=NWRES
 IWRS=NIWRES
 IPRD=NIPRED
 IRS=NIRES
 
 MONITORING OF SEARCH:

 
0ITERATION NO.:    0    OBJECTIVE VALUE:   113.813298695168        NO. OF FUNC. EVALS.:   6
 CUMULATIVE NO. OF FUNC. EVALS.:        6
 NPARAMETR:  1.9600E+00  8.6000E-02  4.0000E-02  1.5190E+00  3.7600E-03 -6.6000E-02  1.3300E-04  7.5200E-03  4.8100E-01  4.7600E-01

 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01 -1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01

 GRADIENT:   7.6043E+00 -5.2852E+01  3.6904E+01  1.7823E+00 -1.3404E+01  5.3011E+00  2.5823E+00 -1.5609E+01 -2.5459E-03 -1.9568E+01

 
0ITERATION NO.:    1    OBJECTIVE VALUE:   112.906821862204        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:       14
 NPARAMETR:  1.9558E+00  8.7300E-02  3.9583E-02  1.5175E+00  3.9011E-03 -6.6959E-02  1.3354E-04  7.8457E-03  5.2431E-01  4.8132E-01

 PARAMETER:  9.7842E-02  1.1500E-01  8.9526E-02  9.9494E-02  1.0380E-01 -1.0150E-01  9.9267E-02  1.0443E-01  1.0000E-01  1.0555E-01

 GRADIENT:   6.5403E+00 -2.2504E+01  1.1281E+01  2.3989E+00 -8.8758E+00  3.5747E+00 -6.2804E-01  3.6952E+01  1.7571E-02 -1.5463E+01

 
0ITERATION NO.:    2    OBJECTIVE VALUE:   112.731466605182        NO. OF FUNC. EVALS.:  10
 CUMULATIVE NO. OF FUNC. EVALS.:       24
 NPARAMETR:  1.9482E+00  8.8738E-02  3.9214E-02  1.5135E+00  4.1065E-03 -6.8353E-02  1.3465E-04  6.6543E-03  3.8273E-01  4.9040E-01

 PARAMETER:  9.3963E-02  1.3134E-01  8.0166E-02  9.8184E-02  1.0941E-01 -1.0375E-01  9.9268E-02  8.9086E-02  9.9993E-02  1.1490E-01

 GRADIENT:   6.8246E+00 -2.1790E+00 -4.7457E+00  3.2425E+00 -1.0401E+01  5.7470E+00  3.7243E+00 -9.0694E+01 -1.3609E-02 -1.1353E+01

 
0ITERATION NO.:    3    OBJECTIVE VALUE:   112.180159468417        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:       32
 NPARAMETR:  1.9088E+00  8.9291E-02  3.9939E-02  1.4841E+00  5.0161E-03 -7.5953E-02  1.4040E-04  6.5313E-03  3.7802E-01  5.3061E-01

 PARAMETER:  7.3526E-02  1.3756E-01  9.8466E-02  8.8363E-02  1.3497E-01 -1.1643E-01  9.8991E-02  8.8436E-02  9.9956E-02  1.5430E-01

 GRADIENT:   5.5995E+00 -8.1992E+00  6.8221E+00  3.6439E+00 -5.0550E+00  3.3207E+00  3.6889E+00 -9.0984E+01 -1.1813E-02  1.8529E+00

 
0ITERATION NO.:    4    OBJECTIVE VALUE:   111.993266021084        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       39
 NPARAMETR:  1.7095E+00  9.0193E-02  4.0424E-02  1.2640E+00  6.1497E-03 -8.9019E-02  1.5460E-04  6.1892E-03  3.5889E-01  4.9241E-01

 PARAMETER: -3.6763E-02  1.4760E-01  1.1054E-01  8.1005E-03  1.7930E-01 -1.4786E-01  1.0396E-01  8.5849E-02  9.9630E-02  1.1695E-01

 GRADIENT:   9.4930E-01 -4.2542E+00  5.1288E+00  2.8230E+00  3.9659E+00 -6.8014E-01  3.6745E+00 -1.1176E+02 -1.1954E-02 -1.2735E+01

 
0ITERATION NO.:    5    OBJECTIVE VALUE:   111.739715080065        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       46
 NPARAMETR:  1.6264E+00  9.0049E-02  4.0292E-02  1.1469E+00  5.2614E-03 -8.4091E-02  1.5045E-04  6.3217E-03  3.6323E-01  4.9916E-01

 PARAMETER: -8.6550E-02  1.4600E-01  1.0728E-01 -4.0502E-02  1.6104E-01 -1.4663E-01  1.1046E-01  8.6390E-02  9.9398E-02  1.2376E-01

 GRADIENT:  -1.7298E+00 -1.3528E+00  6.5199E-01  2.3030E+00  7.0406E-01  5.0577E-01  3.3243E+00 -1.0287E+02 -9.6109E-03 -1.0544E+01

 
0ITERATION NO.:    6    OBJECTIVE VALUE:   111.532469247037        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       53
 NPARAMETR:  1.6605E+00  9.1436E-02  4.0830E-02  8.9945E-01  4.5542E-03 -7.7186E-02  1.5548E-04  6.5807E-03  3.7453E-01  5.0152E-01

 PARAMETER: -6.5798E-02  1.6129E-01  1.2053E-01 -1.6201E-01  1.5740E-01 -1.5198E-01  1.3408E-01  8.7695E-02  9.8620E-02  1.2612E-01

 GRADIENT:   7.0149E-01 -1.7856E+00  6.9708E+00 -4.9551E-01  2.6139E-01  5.1244E-01  3.1451E+00 -9.4317E+01 -6.6951E-03 -1.1030E+01

 
0ITERATION NO.:    7    OBJECTIVE VALUE:   111.509666400812        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:       61
 NPARAMETR:  1.6662E+00  9.0821E-02  4.0512E-02  8.7347E-01  4.6096E-03 -7.7335E-02  1.5809E-04  6.6046E-03  3.7536E-01  5.0407E-01

 PARAMETER: -6.2384E-02  1.5454E-01  1.1271E-01 -1.7667E-01  1.6167E-01 -1.5452E-01  1.3914E-01  8.7768E-02  9.8446E-02  1.2864E-01

 GRADIENT:   8.8859E-01 -3.2298E+00  4.7521E+00 -8.7459E-01  1.2995E+00 -2.2590E-02  3.2271E+00 -9.3469E+01 -6.2946E-03 -1.0639E+01

 
0ITERATION NO.:    8    OBJECTIVE VALUE:   111.323547432968        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       68
 NPARAMETR:  1.6613E+00  9.0985E-02  4.0526E-02  8.9294E-01  4.8452E-03 -7.0721E-02  1.8144E-04  7.5886E-03  4.1613E-01  5.0768E-01

 PARAMETER: -6.5323E-02  1.5634E-01  1.1305E-01 -1.6565E-01  1.6807E-01 -1.3976E-01  2.1328E-01  9.2648E-02  9.4841E-02  1.3222E-01

 GRADIENT:   5.6799E-01 -1.7781E+00  3.1480E+00 -4.2218E-01  1.2879E+00  2.4051E-01  3.3012E+00 -6.7156E+01  1.6781E-03 -8.9853E+00

 
0ITERATION NO.:    9    OBJECTIVE VALUE:   111.065031386817        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       75
 NPARAMETR:  1.6715E+00  9.0561E-02  4.0405E-02  9.2133E-01  4.3634E-03 -8.6440E-02  1.6086E-04  7.2716E-03  4.2977E-01  5.1742E-01

 PARAMETER: -5.9225E-02  1.5167E-01  1.1007E-01 -1.4999E-01  1.4901E-01 -1.6817E-01  1.6262E-01  9.3900E-02  8.7109E-02  1.4171E-01

 GRADIENT:   4.9053E-01 -5.4396E-01  9.3287E-01 -1.9272E-01  5.2990E-01 -1.1327E-01  1.2319E+00 -3.1090E+01  1.0530E-02 -5.1426E+00

 
0ITERATION NO.:   10    OBJECTIVE VALUE:   111.035486991641        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:       82
 NPARAMETR:  1.6613E+00  9.0042E-02  4.0328E-02  9.0583E-01  5.3175E-03 -4.0735E-02  1.6091E-04  7.2110E-03  4.3063E-01  5.2583E-01

 PARAMETER: -6.5350E-02  1.4593E-01  1.0817E-01 -1.5848E-01  1.8314E-01 -7.9924E-02  1.2368E-01  9.4695E-02  7.4119E-02  1.4978E-01

 GRADIENT:   9.8448E-02 -2.2232E+00  1.8714E+00 -2.0012E-01  6.2178E-01  6.6975E-01  5.3977E-01 -1.5876E+01  1.2418E-02 -2.6279E+00

 
0ITERATION NO.:   11    OBJECTIVE VALUE:   110.977727116211        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.6607E+00  9.0251E-02  4.0357E-02  9.2299E-01  4.3551E-03 -7.9908E-02  1.5746E-04  7.4432E-03  4.5440E-01  5.3484E-01

 PARAMETER: -6.5679E-02  1.4825E-01  1.0888E-01 -1.4910E-01  1.4859E-01 -1.5532E-01  1.5077E-01  9.6743E-02  5.5347E-02  1.5827E-01

 GRADIENT:   6.6854E-01 -8.6856E-01  7.1749E-01  3.3622E-02 -3.3014E-01  1.6108E-01  1.5198E-01 -2.8232E+00  1.5044E-02  5.9901E-01

 
0ITERATION NO.:   12    OBJECTIVE VALUE:   110.977393315371        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:       98
 NPARAMETR:  1.6562E+00  9.0274E-02  4.0367E-02  9.2118E-01  4.3383E-03 -8.0818E-02  1.5751E-04  7.4515E-03  4.5536E-01  5.3459E-01

 PARAMETER: -6.8421E-02  1.4850E-01  1.0912E-01 -1.5008E-01  1.4816E-01 -1.5724E-01  1.5139E-01  9.6829E-02  5.1457E-02  1.5804E-01

 GRADIENT:   3.9849E-01 -6.0478E-01  5.1430E-01  3.6533E-02 -2.9902E-01  1.2222E-01  1.1061E-01 -2.0047E+00  1.5158E-02  5.0144E-01

 
0ITERATION NO.:   13    OBJECTIVE VALUE:   110.976546827613        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:      106
 NPARAMETR:  1.6548E+00  9.0283E-02  4.0370E-02  9.2007E-01  4.3468E-03 -8.0755E-02  1.5758E-04  7.4523E-03  4.5566E-01  5.3430E-01

 PARAMETER: -6.9265E-02  1.4860E-01  1.0920E-01 -1.5068E-01  1.4854E-01 -1.5722E-01  1.5126E-01  9.6864E-02  4.7121E-02  1.5777E-01

 GRADIENT:  -3.4264E-01 -4.1319E-01  3.7230E-01  3.1988E-02 -2.3901E-01  9.6322E-02  8.5385E-02 -1.5589E+00  1.5158E-02  3.9523E-01

 
0ITERATION NO.:   14    OBJECTIVE VALUE:   110.976120588918        NO. OF FUNC. EVALS.:   7
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  1.6555E+00  9.0322E-02  4.0373E-02  9.1847E-01  4.3834E-03 -8.0632E-02  1.5762E-04  7.4488E-03  4.5675E-01  5.3292E-01

 PARAMETER: -6.8842E-02  1.4904E-01  1.0928E-01 -1.5155E-01  1.4992E-01 -1.5711E-01  1.4998E-01  9.6985E-02  3.1027E-02  1.5648E-01

 GRADIENT:   4.4973E-01  3.2019E-02  3.6643E-02  2.7637E-04  6.1658E-02 -2.5767E-02 -1.2247E-02  1.6383E-01  1.5148E-02 -5.6686E-02

 
0ITERATION NO.:   15    OBJECTIVE VALUE:   110.976120588918        NO. OF FUNC. EVALS.:  17
 CUMULATIVE NO. OF FUNC. EVALS.:      130
 NPARAMETR:  1.6555E+00  9.0322E-02  4.0373E-02  9.1847E-01  4.3834E-03 -8.0632E-02  1.5762E-04  7.4488E-03  4.5675E-01  5.3292E-01

 PARAMETER: -6.8842E-02  1.4904E-01  1.0928E-01 -1.5155E-01  1.4992E-01 -1.5711E-01  1.4998E-01  9.6985E-02  3.1027E-02  1.5648E-01

 GRADIENT:  -2.0597E-01 -4.9415E-01 -7.1803E-01  2.7637E-04  6.1658E-02 -2.5767E-02 -1.2247E-02  1.6383E-01  1.5148E-02 -9.0499E-02

 
0ITERATION NO.:   16    OBJECTIVE VALUE:   110.974107924916        NO. OF FUNC. EVALS.:  14
 CUMULATIVE NO. OF FUNC. EVALS.:      144
 NPARAMETR:  1.6556E+00  9.0438E-02  4.0424E-02  9.1757E-01  4.3699E-03 -8.1109E-02  1.5792E-04  7.4638E-03  4.5737E-01  5.3292E-01

 PARAMETER: -6.8769E-02  1.5032E-01  1.1055E-01 -1.5204E-01  1.4953E-01 -1.5812E-01  1.5148E-01  9.7043E-02  2.1870E-02  1.5648E-01

 GRADIENT:  -6.6198E-01 -1.2754E-01 -5.2679E-01 -1.1279E-03  4.0037E-02 -1.8460E-02 -8.3056E-03  1.3813E-01  1.4876E-02 -6.9399E-02

 
0ITERATION NO.:   17    OBJECTIVE VALUE:   110.972570737394        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.6569E+00  9.0588E-02  4.0495E-02  9.1762E-01  4.3488E-03 -8.1755E-02  1.5815E-04  7.4805E-03  4.5809E-01  5.3309E-01

 PARAMETER: -6.7972E-02  1.5198E-01  1.1231E-01 -1.5201E-01  1.4881E-01 -1.5938E-01  1.5305E-01  9.7112E-02  3.4135E-05  1.5663E-01

 GRADIENT:  -1.0099E+00  4.2427E-01 -1.1332E-01  5.0120E-04 -1.8318E-02  5.4518E-03 -5.6543E-03  1.1651E-01  1.4243E-02  2.2187E-02

 
0ITERATION NO.:   18    OBJECTIVE VALUE:   110.972570737394        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.6571E+00  9.0585E-02  4.0495E-02  9.1767E-01  4.3478E-03 -8.1751E-02  1.5812E-04  7.4811E-03  4.5819E-01  5.3299E-01

 PARAMETER: -6.7972E-02  1.5198E-01  1.1231E-01 -1.5201E-01  1.4881E-01 -1.5938E-01  1.5305E-01  9.7112E-02  3.4135E-05  1.5663E-01

 GRADIENT:  -2.0935E+00  1.2150E+00  2.8313E-01 -4.3041E-01  1.1295E+00 -1.2681E-01  1.5423E+00 -1.0952E+00 -1.1195E+00  5.4912E-01

 NUMSIGDIG:         2.9         3.5         4.4         3.7         3.6         4.1         3.4         3.9         3.5         3.2

 
 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      181
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:        -1.6234E-02 -2.2335E-04 -7.1104E-03
 SE:             2.6555E-01  3.3141E-03  1.8181E-01
 N:                      12          12          12
 
 P VAL.:         9.5125E-01  9.4627E-01  9.6880E-01
 
 ETASHRINKSD(%)  1.0000E-10  4.6518E+00  2.8106E+00
 ETASHRINKVR(%)  1.0000E-10  9.0873E+00  5.5421E+00
 EBVSHRINKSD(%)  2.8380E+00  7.3200E+00  5.8930E+00
 EBVSHRINKVR(%)  5.5955E+00  1.4104E+01  1.1439E+01
 EPSSHRINKSD(%)  1.0427E+01
 EPSSHRINKVR(%)  1.9767E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          108
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    198.490723172209     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:    110.972570737394     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       309.463293909604     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                            36
  
 #TERE:
 Elapsed estimation  time in seconds:    10.71
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.07
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                        FIRST ORDER CONDITIONAL ESTIMATION                      ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************      110.973       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                        FIRST ORDER CONDITIONAL ESTIMATION                      ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         1.66E+00  9.06E-02  4.05E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.18E-01
 
 ETA2
+        4.35E-03  1.58E-04
 
 ETA3
+       -8.18E-02  7.48E-03  4.58E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        5.33E-01
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.58E-01
 
 ETA2
+        3.61E-01  1.26E-02
 
 ETA3
+       -1.26E-01  8.79E-01  6.77E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        7.30E-01
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                        FIRST ORDER CONDITIONAL ESTIMATION                      ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         1.51E-01  6.42E-03  3.07E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        5.56E-01
 
 ETA2
+        5.71E-03  1.41E-04
 
 ETA3
+        2.28E-01  5.36E-03  2.26E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.31E-02
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        2.90E-01
 
 ETA2
+        4.57E-01  5.59E-03
 
 ETA3
+        3.40E-01  1.44E-01  1.67E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        8.95E-03
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                        FIRST ORDER CONDITIONAL ESTIMATION                      ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11      OM12      OM13      OM22      OM23      OM33      SG11  
 
 TH 1
+        2.27E-02
 
 TH 2
+        3.43E-04  4.12E-05
 
 TH 3
+        1.81E-04  1.48E-05  9.40E-06
 
 OM11
+        1.01E-02 -2.60E-04  1.10E-04  3.09E-01
 
 OM12
+       -2.00E-05 -4.18E-06 -3.52E-07  2.98E-04  3.27E-05
 
 OM13
+       -2.31E-03 -1.80E-04 -7.03E-05 -2.93E-02  1.01E-03  5.22E-02
 
 OM22
+       -4.62E-07  1.49E-07  1.47E-08 -7.85E-06  2.15E-07  4.25E-06  1.98E-08
 
 OM23
+        2.60E-06  6.13E-06  8.16E-07 -1.49E-04 -8.91E-07 -5.14E-05  6.88E-07  2.87E-05
 
 OM33
+        1.32E-03  2.03E-04  6.64E-05  1.70E-03 -2.34E-04 -1.20E-02  2.07E-05  1.05E-03  5.12E-02
 
 SG11
+       -2.47E-05 -1.12E-06 -9.01E-07 -3.64E-04 -6.77E-06 -1.53E-04 -4.08E-08 -1.44E-06  4.72E-06  1.71E-04
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                        FIRST ORDER CONDITIONAL ESTIMATION                      ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11      OM12      OM13      OM22      OM23      OM33      SG11  
 
 TH 1
+        1.51E-01
 
 TH 2
+        3.54E-01  6.42E-03
 
 TH 3
+        3.92E-01  7.50E-01  3.07E-03
 
 OM11
+        1.20E-01 -7.30E-02  6.47E-02  5.56E-01
 
 OM12
+       -2.32E-02 -1.14E-01 -2.01E-02  9.39E-02  5.71E-03
 
 OM13
+       -6.71E-02 -1.23E-01 -1.00E-01 -2.31E-01  7.74E-01  2.28E-01
 
 OM22
+       -2.18E-02  1.65E-01  3.40E-02 -1.00E-01  2.68E-01  1.32E-01  1.41E-04
 
 OM23
+        3.23E-03  1.78E-01  4.97E-02 -4.99E-02 -2.91E-02 -4.20E-02  9.12E-01  5.36E-03
 
 OM33
+        3.87E-02  1.40E-01  9.57E-02  1.35E-02 -1.81E-01 -2.32E-01  6.49E-01  8.66E-01  2.26E-01
 
 SG11
+       -1.25E-02 -1.34E-02 -2.25E-02 -5.01E-02 -9.06E-02 -5.12E-02 -2.22E-02 -2.06E-02  1.59E-03  1.31E-02
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                        FIRST ORDER CONDITIONAL ESTIMATION                      ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11      OM12      OM13      OM22      OM23      OM33      SG11  
 
 TH 1
+        5.37E+01
 
 TH 2
+       -2.44E+02  6.52E+04
 
 TH 3
+       -6.34E+02 -9.89E+04  2.78E+05
 
 OM11
+       -1.81E+00  8.83E+01 -1.14E+02  1.22E+02
 
 OM12
+        3.60E+01  5.09E+03 -1.51E+04 -3.75E+04  1.19E+07
 
 OM13
+       -1.73E+00  1.24E+02  5.94E+01  6.57E+02 -2.08E+05  3.66E+03
 
 OM22
+        1.69E+03  9.35E+04 -3.33E+05  2.96E+06 -9.40E+08  1.64E+07  7.48E+10
 
 OM23
+        1.21E+02 -3.36E+04  6.23E+04 -1.03E+05  3.28E+07 -5.73E+05 -2.61E+09  9.13E+07
 
 OM33
+       -2.94E+00  5.77E+02 -1.15E+03  9.02E+02 -2.86E+05  5.01E+03  2.28E+07 -7.99E+05  7.04E+03
 
 SG11
+        3.33E-01  9.58E+01  4.17E+02 -8.29E+02  2.66E+05 -4.63E+03 -2.11E+07  7.37E+05 -6.44E+03  1.19E+04
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       15.413
Stop Time: 
Fri 05/18/2018 
03:36 PM