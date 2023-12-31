


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        SIZES FSIZES PRSIZES                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: MODULE for NONMEM and its components.
 CONTEXT: Source code

 DISCUSSION:
 File  SIZES is supplied on the NONMEM distribution media and is copied
 to the resources directory.  It contains FORTRAN PARAMETER  statements
 giving  the default values of symbolic parameters used in source code.
 Constants from SIZES may be incorporated into the source code by means
 of  the FORTRAN USE statement.  Some of these constants describe array
 sizes.

 With NONMEM 7.2 and higher, there are  several  changes  vs.  previous
 releases.

      There  is  only  one version of SIZES (i.e., there is no longer a
      SIZES_reg or SIZES_big).

      Many NONMEM and NM-TRAN arrays are allocated dynamically  at  run
      time.

      NM-TRAN  creates  a  subroutine FSIZESR for NONMEM.  FSIZESR con- |
      tains values for some parameters (e.g.,  LVR)  that  are  exactly |
      what  is  needed  for  the  current problem, and contains 0's for |
      parameters that NM-TRAN cannot assess.  Default values for param- |
      eters  that are 0 in FSIZESR are obtained by NONMEM from the file |
      SIZES.f90.                                                        |

      Note that file FSIZES is a  convenient  reference  for  users  to |
      view,  but  is  not used.  SUBROUTINE FSIZESR in FSUBS is what is |
      actually used during the NONMEM run.

      NM-TRAN creates a file prsizes.f90 for PREDPP.  PREDPP arrays are
      allocated  statically  but  may  be re-compiled at run time using
      parameters defined in prsizes.f90.  Some constants in prsizes are
      assessed  for the current problem; others are copied from default
      values in SIZES.f90.

      The user may override many of the parameter values in FSIZES  and
      prsizes  with  the  $SIZES record.  Any value specified by $SIZES
      will override both the default in SIZES.f90 and  the  value  that
      NM-TRAN would have specified.
      (See $sizes)
      As  of NONMEM 7.3, NMTRAN determines the maximum number of obser- |
      vation records (MDV=0) that occur in any subject, among all  data |
      files  used  in the entire control stream file.  If this value is |
      greater than the NO value listed in SIZES.f90, it will set NO  to |
      this larger size.  Thus, users no longer have to be conscientious |
      of sizing the NO parameter. However, there is no  guarantee  that |
      NMTRAN  will correctly assess NO for the entire scope of the con- |
      trol stream file for all types of problems.  Should  this  occur, |
      NONMEM  may  issue an error, and the user will need to set the NO |
      value with a $SIZES record.

 Constants that may be changed with $SIZES record:

 LTH LVR LVR2 LPAR MMX MAXIDS NO PD MAXOMEG MAXPTHETA

 Consants that are also in prsizes:

 PC PCT PIR PD PAL MAXFCN

 Constants that may not be changed using $SIZES record

 LSTEXT MAXXNAME MAXVRESWRES MAXVRESN MXNAME DIMQ PPR  PW  SD  SCO  SDF
 MAX_EXTRA

 The  following  is the descriptive comment and value of each parameter
 in the version of SIZES that is supplied on the NONMEM 7  distribution
 medium.  More information including descriptions of the buffers may be
 found in the INTRODUCTION TO NONMEM 7.

 NLUSER=100
      Maximum number of logical I/O units user may use.

 LTH=100
      MAX. NO. OF THETA'S.  Dynamically  sized,  or  set  by  user  via
      $SIZES record

 LVR=30
      MAX. NO. OF ETA'S + EPS'S.  Dynamically sized, or set by user via
      $SIZES record

 LVR2=20
      MAX. NO. OF ETA'S PERMITTED WHEN LAPLACIAN METHOD IS USED.  Value
      may  be  over-ridden by user via $SIZES record Value may be over-
      ridden by user via $SIZES record

 NO=250
      MAX NO. OF OBSERVATION RECORDS / INDIVIDUAL RECORD

 MMX=10
      MAX NO. OF MIXTURE SUBPOPULATIONS.  Dynamically sized, or set  by
      user via $SIZES record

 LNP4=4000
      SIZE  OF  COMMON  NMPRD4.   Value  may be over-ridden by user via
      $SIZES record

 LSUPP=4050
      MAX. NO. OF POINTS OF SUPPORT WITH NONPARAMETRIC ESTIMATE.  Value
      may be over-ridden by user via $SIZES record

 LIM7=2
      SIZE  OF  BUFFER  7.  DO NOT GO LOWER THAN 2.  Value may be over-
      ridden by user via $SIZES record

 LWS3=9000
      SIZE OF WORKING SPACE 3 AT LEAST AS LARGE AS: NS*NOETAS**2, WHERE
      NS  IS  THE NO. OF DIRECTIONS USED WITH THE STIELTJES METHOD, AND
      NOETAS IS THE NUMBER OF ETA'S.  Value may be over-ridden by  user
      via $SIZES record

 MAXIDS=10000
      MAX.  NO.  OF INDIVIDUALS IN DATA SET.  Dynamically sized, or set
      by user via $SIZES record

 LIM1=10000
      SIZE OF BUFFER 1.  Value may be over-ridden by  user  via  $SIZES
      record

 LIM2=100000
      SIZE  OF  BUFFER  2.  Value may be over-ridden by user via $SIZES
      record

 LIM3=10000
      SIZE OF BUFFER 3. DO NOT GO LOWER THAN 2.  Value may be over-rid-
      den by user via $SIZES record

 LIM4=1000
      SIZE  OF  BUFFER  4: LIM4=NUMBER OF SUBJECTS.  Value may be over-
      ridden by user via $SIZES record

 LIM5=200
      SIZE OF BUFFER 5.  Value may be over-ridden by  user  via  $SIZES
      record

 LIM6=400
      SIZE  OF  BUFFER  6.  Value may be over-ridden by user via $SIZES
      record

 LIM8=200
      SIZE OF BUFFER 8.  Value may be over-ridden by  user  via  $SIZES
      record

 LIM10=100000
      SIZE  OF  BUFFER 10.  Value may be over-ridden by user via $SIZES
      record

 LIM11=25
      SIZE OF BUFFER 11.  Value may be over-ridden by user  via  $SIZES
      record

 LIM13=1000
      SIZE  OF BUFFER 13: LIM13=NUMBER OF SUBJECTS.  Value may be over-
      ridden by user via $SIZES record

 LIM15=1000
      SIZE OF BUFFER 15: LIM15=NUMBER OF SUBJECTS.  Value may be  over-
      ridden by user via $SIZES record

 LIM16=400
      SIZE  OF  BUFFER 16.  Value may be over-ridden by user via $SIZES
      record

 LSTEXT=70000
      AT LEAST (MAXOMEG*(MAXOMEG+1)/2+MAXPHETA+5)*NUMTEXT.   LSTEXT  is
      maximun  number  of characters to a single line of the raw output
      file specified by $EST FILE=.  where NUMTEXT is number of charac-
      ters  needed  to represent a number and its delimiter.  For exam-
      ple, FORMAT=,1PE12.5 takes up NUMTEXT=13 characters.

 LSFORM=2048
      LSFORM is the character length of FORMAT TFORMATL,RFORMATL,  per-
      taining  to  the  full  length of LFORMAT, RFORMAT for the $TABLE
      record.

 MAXRECID=200
      Maximum number of records in any individual record.  For  use  by
      PREDPP.  Value may be over-ridden by user via $SIZES record

 PC=30
      MAX. NO. OF COMPARTMENTS MAXIMUM IS 99.  Value may be over-ridden
      by user via $SIZES record

 PCT=30
      MAX. NO. OF MODEL EVENT TIMES.  Value may be over-ridden by  user
      via $SIZES record

 MAXNRDS=PC
      MAXNRDS: MAX. NO. OF DELAY COMPARTMENTS FOR ADVAN16, ADVAN17, and
      ADVAN18.  YOU MAY WANT TO MAKE IT LOWER THAN PC AND  SAVE  MEMORY
      Value may be over-ridden by user via $SIZES record

 PAST_SIZE=4000
      RESOLUTION  (NUMBER  OF  DETAILED  POINTS) FOR DELAY COMPARTMENTS
      STORAGE FOR ADVAN16 AND ADVAN17.  Value  may  be  over-ridden  by
      user via $SIZES record

 PCT_BIG=10000
      MAX.  NO.  OF  MODEL EVENT TIMES THAT NMTRAN CAN PROCESS.^M ! PG:
      SIZE OF GG; MAX. NO. OF BASIC+ADDITIONAL PK PARAMS.^M !     (MAX-
      IMUM IS PCT+99)^M

 PIR=700
      SIZE  OF  COMPACT  DA/DP/DT  ARRAYS.  Value may be over-ridden by
      user via $SIZES record

 PD=50
      CHANGED TO INCREASE DATA ITEMS PER DATA RECORD FROM 20 TO  50  IT
      IS ALSO SIZE OF VDATREC DATA ARRAY. THIS IS DONE TO SEPARATE VDA-
      TREC AND VRESWRES VARIABLES/ Dynamically sized, or  set  by  user
      via $SIZES record.

 MAXIC=90
      MAXIC:  MAXIMUM NUMBER OF ACTIVE INFUSIONS FOR PREDPP.  Value may
      be over-ridden by user via $SIZES record

 PDT=500
      MAXIMUM NUMBER OF TABLE ITEMS/PRED-DEFINED ITEMS.  Value  may  be
      over-ridden by user via $SIZES record

 PAL=50
      NO.  OF ADDITIONAL AND LAGGED DOSES.  Value may be over-ridden by
      user via $SIZES record

 MAXFCN=1000000
      MAX. NO. OF CALLS IN GENERAL NON-LINEAR  MODELS  IMAX  IN  MODULE
      PRCOM_INT OVERRIDES.  Value may be over-ridden by user via $SIZES
      record

 STIELTJ_SIZE=101

      The next four are for internal use, pertaining to  various  addi-
      tional weighted residual diagnostics

 MAXXNAME=40

 MAXVRESWRES=39

 MAXVRESN=9

 MXNAME=40

 DIMTMP=500
      RELATED  TO  THE  NUMBER OF USER-DEFINED VARIABLES.  Value may be
      over-ridden by user via $SIZES record

 DIMCNS=500
      RELATED TO THE TOTAL NUMBER OF CONSTANTS.  Value may be over-rid-
      den by user via $SIZES record

 DIMNEW=1000
      RELATED TO THE TOTAL NUMBER OF INTERMEDIATE VARIABLES.  Value may
      be over-ridden by user via $SIZES record

 DIMQ=99999
      ARRAY SIZE FOR LOGICAL CONDITIONS

 FL=49
      LOGICAL UNIT NUMBER FOR FLIB

 DIMVRB=1000
      MAX. NO. OF LINES OF VERBATIM CODE

 PL=10
      MAXIMUM DEPTH OF NESTED IF STATEMENTS

 NFUNCX=100
      MAXIMUM NUMBER OF USER FUNCTIONS.

 NVECX=100
      MAXIMUM NUMBER OF USER VECTORS.

 MAXOTHER=1000
      Maximum number of filenames listed on $SUBROUTINE OTHER=filename

 SD=30
      LENGTH OF DATA LABEL

 FSD=67000
      LENGTH OF CONTROL STREAM FILE STRING (ALSO, MULTIPLE  LINES  CON-
      CATENATED WITH & MAY NOT EXCEED FSD)

 FSD1=67001
      LENGTH OF CONTROL STREAM FILE STRING

 SCO=30
      STRING LENGTH OF NUMBER IN $THETA, $OMEGA, $SGIGMA RECORDS.

 SDF=24
      LENGTH OF DATA ITEM

 NFSIZES=50
      NUMBER OF ITEMS LISTED IN FSIZES FILE

 NPRSIZES=14
      NUMBER OF PREDPP ITEMS IN PRSIZES.F90

 CONSTANTS FOR Monte Carlo, EM methods

 MAX_EXTRA=20
      Numer of $EST statements allowed per problem

 NPOPMIXMAX=10
      Now  dynamically  sized  to MAXIMUM NUMBER OF SUB-POPULATIONS FOR
      MIXTURE MODELS MMX.  Or, set by user via $SIZES record

 MAXOMEG=70
      Now dynamically sizes to LVR, OR MAXIMUM  OMEGA  DIMENSION  ETAS.
      Or, set by user via $SIZES record (but should be left alone)

 MAXPTHETA=90
      Now dynamically sized to LTH, MAXIMUM NUMBER OF THETA, PLUS LOWER
      TRIANGLE OF SIGMA.  THUS, IF NUMBER OF THETAS IS N, AND DIMENSION
      OF  SIGMAS  IS  M, THEN NEED MAXPTHETA=N + M*(M+1)/2.  OR, may be
      set by user via $SIZES record (but should be left alone).

 MAXITER=210
      Maximum number of previous iterations to incorporate  into  Monte
      Carlo convergence tests Effective CITER is <=MAXITER

 ISAMPLEMAX=10
      FOR SAEM METHOD, EFFECTIVE ISAMPLE<=ISAMPLEMAX

 MAXSIDS=100
      MAXIMUM NUMER OF SUPER ID ITEMS IN DATA SET

 MAXSIDL=0
      MAXIMUM  LEVEL  OF  SUPER  IDS.   MAXSIDL=0  MEANS  NO SUPER IDS.
      Default 0, until $LEVEL is used, in which case it will be dynami-
      cally sized.

 MAXFTEXT=100
      Maximum of PRDERR message lines

 PNM_MAXNODES=100
      MAXIMUM NUMBER OF PARALLELIZATION NODES

 PNM_BUFFER_SIZE=100000
      INTERNAL BUFFER SIZE, FOR EFFICIENT PACKAGING AND SENDING BETWEEN
      PROCESSES

 REFERENCES: Guide III, section V.3.0 
 REFERENCES: Guide Introduction_7
