


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            $ABBREVIATED                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Provides instructions about abbreviated code
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $ABBREVIATED  [COMRES=n1] [COMSAV=n2]
               [DERIV2=NO] [DERIV2=NOCOMMON] [DERIV1=NO]
               [FASTDER | NOFASTDER]
               [CHECKMU | NOCHECKMU]
               [DES=COMPACT|DES=FULL]
               [REPLACE left_string = right_string ] ...
               [DECLARE [type] [DOWHILE] name [(dimension [,dimension])] ...
               [PROTECT]
               [FUNCTION function_name(input_vector_name,dimension[,usage])]
               [VECTOR input_vector_name(dimension)]

 SAMPLE:
 $ABBREVIATED    COMRES=2

 DISCUSSION:
 Optional.   May be used when $PK, $ERROR, or $PRED abbreviated code is
 present. Must precede all blocks of  abbreviated  code.   With  NONMEM
 7.4, may also be used when there is no abbreviated code.  For example,
 $ABBR REPLACE may be used for  label  substitution  in  NONMEM  report
 files.

 OPTIONS:

 COMRES=n1 ('common reserve')
      COMRES gives instructions to NM-TRAN.
      Values of n1:

   -1 Do not store any variables in the global variable array, VRBL, in
      the NMPRD4 module.

    0 Store  variables  in  NMPRD4  with  no  reserved  positions  (the
      default)

   n1 Store variables in NMPRD4, but reserve the first n1 positions

      With  abbreviated  code, the Ith position in NMPRD4 is referenced
      by COM(I).

      This option is intended for advanced users of NONMEM, e.g.,  when
      abbreviated  code  is  combined with user-supplied subroutines or
      verbatim code.  A user-supplied subroutine may reserve the  first
      n1 positions in NMPRD4 for its use, in which case the option COM-
      RES should be set to n1 to instruct NM-TRAN to skip  these  posi-
      tions;  the  first position used by NM-TRAN for storing variables
      defined in abbreviated code will be position {n1}+1.

      $TABLE  and  $SCATTER  may  explicitly  reference  variables   in
      reserved  positions  1  through  n1  by  COM(1)  through COM(n1),
      respectively, in addition to listing variables defined in  abbre-
      viated code by name.

      An  individual  block  of abbreviated code (e.g. $PK) may include
      the  pseudo-statement  COMRES=-1,  which  prevents  any  variable
      defined in that particular block from being stored in NMPRD4.

 COMSAV= n2 ('common save')
      Values  of  variables  displayed  in  tables and scatterplots are
      obtained from NMPRD4.   There  are  particular  times  when  data
      records  are  passed  to  PRED for the purpose of obtaining these
      values; these are called copying  passes.   The  SAVE  region  of
      NMPRD4  is  the  initial  part  of  this array.  If a variable is
      stored in the SAVE region, then the value of  the  variable  com-
      puted  with  a  given  data  record during a copying pass will be
      found in NMPRD4 when the same record is passed  during  the  next
      copying  pass,  i.e.  it  will  have been saved from the previous
      copying pass.  This is in contrast to the usual behaviour,  where
      with  a  given data record, the value in NMPRD4 is the value com-
      puted with the previous data record.
      n2 is the initial size of the SAVE region,  i.e.  the  number  of
      positions  in  this  region.  n2 =0 is the default value.  n2 may
      not exceed n1.
      The SAVE region has size n2 initially, but NM-TRAN may extend  it
      if  SAVE variables are used.  However, if n2 =-1, the SAVE region
      is not to be extended, and there is to be no  SAVE  region  alto-
      gether.
      (See copying block).

      When  PREDPP is used, and a $PK block is present, NM-TRAN inserts
      code into the PK routine that  stores  the  value  of  COMSAV  at
      ICALL=1.   If  no  $PK  block  is  present, and a $ERROR block is
      present, the code is  inserted  into  the  ERROR  routine.   When
      PREDPP  is  not used, and a $PRED block is present, the generated
      or library PRED routine stores the value of COMSAV at ICALL<=1.

 DERIV1=NO (NM74)
      Prevents the computation of first derivatives.  For example,  may
      be  used when only SAEM or BAYES is performed, or IMPMAP/ITS/FOCE
      are performed using OPTMAP>0 and ETADER>0.  The  global  variable
      NOFIRSTDERCODE is set to 1.

 DERIV2=NO
      Prevents  the computation of second derivatives, which are needed
      only for the Laplacian method.

 DERIV2=NOCOMMON
      Permits the computation of these derivatives, but  prevents  them
      from being stored in the global variable NMPRD4.
      $ESTIMATION METHOD=COND LAPLACIAN may be specified, but variables
      representing second derivatives are not stored in NMPRD4.  There-
      fore,  they  cannot  be displayed in tables and scatterplots.  In
      addition, no variables computed in the $PK block  may  be  refer-
      enced  in  the  $ERROR  block.  This is true whether or not these
      variables happen to have second derivatives, and whether  or  not
      the Laplacian method is used.

 FASTDER, NOFASTDER
      With NONMEM 7.2 and higher, NM-TRAN collects statements that com-
      pute first-partial eta-derivatives together in  FSUBS,  and  they
      are  performed  only when NONMEM sets IFIRSTEM=1.  NOFASTDER pre-
      vents NM-TRAN from doing this, and restores the order  of  state-
      ments in FSUBS to what it was in previous versions.

      FASTDER  requests  that  the  statements be collected, and is the
      default.
      (See Partial Derivative Indicators).

 CHECKMU, NOCHECKMU
      With NONMEM 7.2 and higher, NM-TRAN checks the  MU  model  state- |
      ments  in  abbreviated  code  and issues warning messages if they |
      appear to contain mistakes.  This can take a long time for  large |
      control streams.  Also, in the examples directory, there are con- |
      trol streams for which the check  is  too  difficult  for  NMTRAN |
      (tdist6_sim.ctl  and tdist7.ctl), and some for which the warnings |
      are  inappropriate  (the  "superid"  control  streams  that   use |
      $LEVEL).   NOCHECKMU can be used to prevent NM-TRAN from attempt- |
      ing to check the MU model statements.                             |

      CHECKMU requests that MU model statements be checked, and is  the |
      default.  Neither option affects the generated code.

 DES=FULL
      Arrays of the DES routine are stored in non-compact form.
      With  $ESTIMATION  METHOD=COND LAPLACIAN, the option NUMERICAL is
      also required.
      DES=FULL is the default with  ADVAN9  and  ADVAN15  and  ADVAN17.
      (Prior to NONMEM 7.4, FULL was required with ADVAN13.)

 DES=COMPACT
      Arrays of the DES routine are stored in compact form.
      Required  with Laplacian method; optional otherwise.  This is the
      default, except with ADVAN9 and ADVAN15 and ADVAN17.

 REPLACE left_string = right_string                                     |
      There are several different forms of the REPLACE option.  They do |
      not affect verbatim code.  Case is ignored.                       |

    (1)  Simple replacement                                             |
      REPLACE left_string = right_string                                |

      May  be  used  in  all  blocks of abbreviated code and $TABLE and |
      $SCATTER records.  Left_string is replaced by right_string.   The
      search   is   "anchored"  by  a  Fortran  identifier.   That  is,
      left_string is parsed for the first Fortran  identifier  and  the
      abbreviated  code is searched for the same identifier.  Following
      the identifier, the entire left_string must match the abbreviated
      code.   Single  or  double  quotes  around a string are optional.
      They are necessary if the string contains punctuation such  as  a
      comma.   A  given  line  of  abbreviated  code  may have multiple
      replacements; they are applied in the order that they  appear  in
      the NM-TRAN control stream.

      Examples:

      $ABBR REPLACE PI=3.14159265
      $ABBR REPLACE THETA(CL)=THETA(4)
      $ABBR REPLACE ETA(CL)=ETA(5)
      $ABBR REPLACE K34="3,4"
       ...
      CL=PI*THETA(CL)*EXP(ETA(CL))
      PRINT *,OMEGA(K34)

      The post-replacement code is
      CL=3.14159265*THETA(4)*EXP(ETA(5))
      PRINT *,OMEGA(3,4)

      With  this feature, subscripts of THETA, ETA, and EPS and ERR may |
      be given symbolic names in abbreviated code, and $ABBR REPLACE is |
      used  to  replace them with integer subscripts.  With NONMEM 7.4, |
      this feature also requests label substitution  in  NONMEM  report |
      files and table files.  For example, the following will cause all |
      appearances of "ETA(CL)" to be replaced by "ETA(3)" in the gener- |
      ated  subroutine,  and  will  cause  "ETA(3)"  to  be replaced by |
      "ETA(CL)" in the NONMEM  report  and  any  table  file  in  which |
      "ETA(3)" is  listed.

      $ABBR REPLACE ETA(CL)=ETA(3)                                      |

      Label  substitution  can be turned off for an entire problem with |
      NOSUB=1 option of the $DEFAULT record.  Label substitution can be |
      turned  off for a specific NONMEM task (and all subsequent tasks) |
      with the NOSUB=1 option of $TABLE, $SCAT, and $ESTIMATION.

      Substitutions will never be made in the additional  output  files |
      *.ext,  ,phi,  etc., to maintain their third party software read- |
      ability.

      With this feature, compartment names may be used instead of  com-
      partment numbers.  For example,
      $ABBR REPLACE A(DEPOT)=A(1)
      $ABBR REPLACE DADT(DEPOT)=DADT(1)
        ...
      $DES
       DADT(DEPOT)=-KA*A(DEPOT)

      This  is  called  "explicit" compartment name substitution.  With
      NONMEM 7.5, compare the "implicit" compartment name  substitution
      feature of the $MODEL record.

      See INTRODUCTION  TO  NONMEM  7,  Symbolic Label Substitutions of
      Model Compartments

    (2) Replacement with selection by data item
      $ABBR REPLACE VAR(d)=VAR(n1,n2,..,nk)
      May be used in $PK, $ERROR, $PRED blocks only.  VAR must  one  of
      ETA,  EPS, THETA. (ERR is not permitted.)  d must be a data item.
      The integer value of d (i.e., INT(d)) is used to  select  one  of
      n1,n2,..,nk.   If  any  of the ni is 0, that position is skipped.
      The effective code is:
      IF (d.eq.1) VAR(d)=VAR(n1)
      IF (d.eq.2) VAR(d)=VAR(n2)
       ...
      The actual generated code uses Q-type  indicator  variables,  and
      variables such as THETA_OCC_ in place of THETA(OCC).

      Example: Suppose OCC is a data item that takes values 1 and 2.
      $ABBR REPLACE THETA(OCC)=THETA(4,7)
        ...
      $PK
      TVCL=THETA(OCC)
      The effective code after replacement is
      IF (OCC==1) TVCL=THETA(4)
      IF (OCC==2) TVCL=THETA(7)

      A  user-defined variable can be given a default value in case the
      implied selection is not satisfied.  E.g., suppose  there  are  3
      doses per subject (DOSN=1,2,3) and the model is:
      $PK ...
         IF(DOSN.EQ.1) F1=1
         IF(DOSN.EQ.2) F1=1*EXP(ETA(2))
         IF(DOSN.EQ.3) F1=1*EXP(ETA(3))

      The model can be implemented using implied selection as follows:
      $ABBR REPLACE ETA(DOSN)=ETA(0,2,3)
      $PK ....
        F1=1
        IF (DOSN>1) F1=1*EXP(ETA(DOSN))

    (3) Replacement with selection by data item and parameter
      $ABBR REPLACE VAR(p_d)=VAR(n1,n2,..,nk)
      p  is a user-defined variable in the abbreviated code.  The order
      may be p_d or d_p This allows a given data item d to be used as a
      selection  variable for more than one parameter, with a different
      choice of elements of VAR.

      Example 1:
      $ABBR REPLACE THETA(SID_KA)=THETA(4,6)
      $ABBR REPLACE THETA(SID_CL)=THETA(5,7)
         ...
      $PK
      KA=THETA(SID_KA)
      CL=THETA(SID_CL)
      which is equivalent to
      $PK
      IF (SID==1) KA=THETA(4)
      IF (SID==2) KA=THETA(6)
      IF (SID==1) CL=THETA(5)
      IF (SID==2) CL=THETA(7)

      Example 2:
      $ABBR REPLACE ETA(OCC_CL)=ETA(5,3)
      $ABBR REPLACE ETA(OCC_V) =ETA(6,4)
         ...
      $PK
      CL=TVCL*EXP(ETA(1)+ETA(OCC_CL))
      V =TVV *EXP(ETA(2)+ETA(OCC_V))
      which is equivalent to
      $PK
      IF (OCC==1) CL=TVCL(EXP(ETA(1)+ETA(5))
      IF (OCC==2) CL=TVCL(EXP(ETA(1)+ETA(3))
      IF (OCC==1) V=TVV(EXP(ETA(2)+ETA(6))
      IF (OCC==2) V=TVV(EXP(ETA(1)+ETA(4))

      The number of values specified for the selection data  item  must
      be consistent for all parameters in which it is used.             |

    (4) Replacement for multiple variables                              |
      $ABBR REPLACE VAR(p1,p2,...,pk)=VAR(n1,n2,...,nk)
      This  form  is  new  to NONMEM 7.4.  The pi are character strings
      separated by commas.  A  character  string  may  not  contain  an
      embedded  space.   The lists on both left and right sides must be
      of the same length.  VAR(pi) is replaced by VAR(ni).   For  exam-
      ple,

      $ABBR REPLACE THETA(CL,V1,Q,V2)=THETA(1,2,3,4)

      This is equivalent to:
      $ABBR REPLACE THETA(CL)=THETA(1)
      $ABBR REPLACE THETA(V1)=THETA(2)
      $ABBR REPLACE THETA(Q)=THETA(3)
      $ABBR REPLACE THETA(V2)=THETA(4)

      Label  substition  occurs with this form, as with Simple replace-
      ment (1).

    Short-hand notation
      A short-hand notation may be used to describe a series of  values
      of ni.  An ascending sequence of ni can be described as
       ,start TO end [BY interval]

      BY  is  optional  and must be > 0. Default is 1.  TO is required.
      If end < start, the sequence is ignored.  At least one comma must
      appear, so NMTRAN knows it is a number list, not a variable name. |
      With NONMEM 7.4, the character : may be used instead of  TO,  and |
      the value of BY may be negative.

      EXAMPLES:
      $ABBR REPLACE THETA(SID_KA)=THETA(10:4 by 3) ; order: 10,7,4
      $ABBR REPLACE THETA(SID_KA)=THETA(4 to 10 by -3) ; order: 10,7,4
      $ABBR REPLACE THETA(SID_KA)=THETA(10 to 4) ; order: 10,9,8,7,6,5,4

      $ABBR REPLACE THETA(SID_KA)=THETA(,4 to 13 by 3,25 to 37 by 4)
      is identical to
      $ABBR REPLACE THETA(SID_KA)=THETA(4,7,10,13,25,29,33,37)

    Files FORIG and FREPL
      When $ABBR REPLACE is coded, NM-TRAN produces two files:
      FORIG Contains the original (pre-replacement) abbreviated code,
      and task specification records $TABLE and $SCATTER.               |
      FREPL Contains the new (post-replacement) abbreviated code
      and task specification records.
      These may be helpful for debugging.

 DECLARE [INTEGER] [DOWHILE] name [(dimension [,dimension])] ...
      One  or  names  may  be  coded. They are referred to as "declared
      variables."  IF INTEGER or DOWHILE is  coded,  the  type  of  the
      variable is integer.  Otherwise, the type of the variable is dou-
      ble precision.  If one or two dimensions are declared, the  vari-
      able  being declared is an array.  Each dimension may be a number
      or a variable or  an  expression.   Constants  defined  in  SIZES
      (e.g.,  NO,  LVR)  may  be used.  Multiple DECLARE options may be
      coded. The characters "DECLARE" are  optional  after  the  first.
      Commas  are  ignored, and type and dimensions must be respecified
      as needed.  No other  options  of  $ABBR  may  appear  after  the
      DECLARE  option(s).   Declared  variables  are  global, i.e., are
      defined in all blocks of abbreviated code.  The  number  of  sub-
      scripts  must agree with the number of dimensions in the declara-
      tion.  Declared variables that are not INTEGER or DOWHILE will be
      random variables if they are assigned in a statement whose right-
      side involves ETA's or EPS's.

      EXAMPLE:

      $ABBR DECLARE A,B(10),C(1,NO-1),INTEGER I J

      Only I is INTEGER.

      Variables may be declared as INTEGER or DOWHILE for use  as  sub-
      scripts  of declared arrays or reserved variables that are arrays
      (but not of ETA, ERR, or EPS).

      Variables used as looping indices in DOWHILE statements  must  be
      declared as DOWHILE variables.

      Declared variables are automatically initialized to 0.

      Elements  of  a  declared  array  may be displayed in WRITE/PRINT
      statement, but not the entire array.  E.g., the following is per-
      mitted
      PRINT *,B(1),B(2)
      but not
      PRINT *,B

      Neither individual elements nor the entire array may be listed in
      $TABLE.  The workaround is to code, e.g.,
           ... code assigning values to B(1) etc. ...
         B1=B(1)
       $TABLE B1
      (See dowhile block).

 PROTECT
      With NONMEM 7.4, a series of routines are available that  protect
      against  domain  violations,  divide  by zero, and floating point
      overflows.  Each of these routines start with the letter P,  fol-
      lowed  by the name of the mathematical operation they are to per-
      form.  For example, PLOG is the protective code routine that per-
      forms  the  LOG operation.  With $ABBR PROTECT, NMTRAN will auto-
      matically replace all relevant function names with the P name.

      (See protect functions).

 FUNCTION function_name(input_vector_name,dimension[,usage])
      With all versions of NONMEM since NONMEM VI, user-supplied  func-
      tions  are permitted in abbreviated code with reserved names such
      as FUNCA, FUNCB, ...,  and argument vectors with  reserved  names
      such as VECTRA, VECTRB, .....  With NONMEM 7.4, extended reserved
      names are recognized.  These are FUNCxy and FUNCxyz,  where  each
      of  x,  y, z stands for an alphabetic character A-Z, e.g., FUNCAB
      or FUNCABC.  Similar extended reserved names for vectors are also
      recognized:  e.g, VECTRAB or VECTRABC.  Reserved argument vectors
      and arrays have dimension 9 and each reserved function may appear
      in abbreviated code at most 9 times.
      (See Abbreviated Function).

      In NONMEM 7.4 the $ABBR FUNCTION option allows user-defined func-
      tion names and user-defined argument vector  names.   The  dimen-
      sions  of  the  argument vector and the maximum number of times a
      given function name may appear in abbreviated code is user-speci-
      fied.

      A user-defined function may be declared as follows:
      $ABBR FUNCTION function_name(input_vector_name,dimension,usage)

 function_name
      is the name of the function.

 input_vector_name
      is the name of an input vector that may be used to pass arguments
      to the function.

 dimension
      specifies how many input arguments function_name  will  use,  and
      defines  input_vector_name as a vector with this length.  "Dimen-
      sion" is a property of both the function and of the input vector.

 usage
      is the maximum number of times the function  may  appear  in  the
      abbreviated  code,  that  is, the maximum number of occurances of
      function_name. It is not an error if there are fewer  occurances.
      If  usage  is  omitted,  NMTRAN  will supply the exact number. If
      usage is coded, NMTRAN will generate an error  message  if  func-
      tion_name appears in abbreviated code more than "usage" number of
      times.

      A vector and its length may be declared independently of a  func-
      tion,
      $ABBR VECTOR input_vector_name (dimension)

 input_vector_name
      is the name of an input vector.

 dimension
      The length of the vector.  The dimension of a vector should be no
      less than the dimension of all the functions which  which  it  is
      used.

      Example:

      $ABBR FUNCTION BIVARIATE(VBI,5,3)

      A  vector VBI is defined of length 5.  There is a function called
      BIVARIATE.  When BIVARIATE is used, the value 5 is passed  to  it
      as  argument NDIM.  BIVARIATE uses 5 elements from the input vec-
      tor.  Function BIVARIATE may appear in abbreviated code at most 3
      times.   BIVARIATE should be present in a source code file listed
      in the $SUBROUTINES record:

      $SUBROUTINES OTHER=filename

      Instructions for coding both reserved functions and  user_defined
      functions are in Abbreviated Function.

      Here  is an example of abbreviated code that uses BVI and BIVARI-
      ATE:
      $PK ...
      VBI(1)=RHO
      VBI(2)=5
      VBI(3)=6
      VBI(4)=1 ;***0 = Upper tail as in Drezner & Wesolowsky
               ; 1 = Bottom tail***;
      VBI(5)=1 ;***0 = 3 pt approximation
                ; 1 = 5 point approximation***;
      BV=BIVARIATE(VBI)

      There is a fully worked out example.
      (See bivariate function).
      Files are  in  the  examples  directory:  bivariate.ctl,  bivari-
      ate.csv, bivariate.f90.

      By  default,  there may be up to 100 different user-defined func-
      tions, which includes functions with  reserved  names  and  user-
      defined  names.  There may also be up to 100 different vectors of
      input arguments, which includes vectors with reserved  names  and
      user-defined  names.   Constants NFUNCX and NVECX in SIZES may be
      used to change these numbers.

      It is preferable that a unique vector should be  associated  each
      function to assure that each vector-function pair are set up with
      comparable dimensions.  This is not necessary.  Any  user-defined
      vector  or  function  may  be  used  with  any reserved vector or
      reserved function, subject to the same restriction. For example:

      $ABBR FUNCTION BIVARIATE(VBI,5)
      $ABBR FUNCTION BIVARIATEQ(VQI,10)

      In abbreviated code, any of the following are permitted:
      BVI=BIVARIATE(VBI)
      BVQ=BIVARIATE(VQI)
      BVQI=BIVARIATEQ(VBI)
      BVQQ=BIVARIATEQ(VQI)

      NM-TRAN gives a warning:
      (WARNING  133) DIMENSION OF VECTOR IS LESS THAN
                     WHAT IS EXPECTED BY FUNCTION
      BVQI=BIVARIATEQ(VBI(001),FNC002_1(1,001),FNC002_2(1,1,001),10)

      This is because vector VBI has dimension 5 and  function  BIVARI-
      ATEQ  was declared with dimension 10.  Note that the results will
      be correct if BIVARIATEQ does not use more than  5  positions  in
      the argument vector.

 A  vector  and its length may be declared independently of a function,
 and vice versa.  The asterisk is used as a placeholder, e.g.,

      $ABBR FUNCTION BIVARIATE(*,5)
      ; when BIVARIATE is called, NDIM will be 5
      $ABBR FUNCTION BIVARIATEQ(*,10)
      ; when BIVARIATEQ is called, NDIM will be 10
      $ABBR VECTOR VQI(15)

 The $ABBREV FUNCTION record may be used to override the  default  set-
 tings for any of the reserved function or vectors.  For example,

      $ABBR FUNCTION FUNCA(VECTRA,25,5)

      VECTRA  will  be  defined with length 25, not 9 as for a reserved
      function.  The code for the function should have NDIM as an argu-
      ment:
      FUNCTION VECTRA(X,X1,X2,NDIM)

      NDIM should be used instead of 9 for vector and array dimensions,
      as shown in  the Abbreviated Function help item.

 REFERENCES: Guide IV, section III.B.7 , IV 
 REFERENCES: Guide Introduction_7
