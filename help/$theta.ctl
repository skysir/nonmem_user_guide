


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $THETA                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives initial estimates and bounds for thetas
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $THETA  value1  [value2]  [value3] ...
         [(value)xn] ...
         [label=value] ...
         [NAMES (label ...)value ...]
         [NUMBERPOINTS=n]
         [ABORT|NOABORT|NOABORTFIRST]

 SAMPLE:
 $THETA  (0,3) 2 FIXED (0,.6,1) 10 (-INF,-2.7,0) (37 FIXED)

 DISCUSSION:
 Gives  initial  estimates and bounds for elements of the THETA matrix.
 Thetas are numbered in the order in which they are defined.

 OPTIONS:

 Each value defines a theta and gives its initial estimate and  bounds.
 A value has one of 4 forms:

 init [FIXED]
      Init is the initial estimate.  If FIXED is used, the final param-
      eter estimate is to be constrained to equal the initial parameter
      estimate.

 ([low,] init [,up] [FIXED])
      Low and up are lower and upper bounds respectively.  They are the
      boundaries for the minimization search.  Commas are optional.  If
      an  upper  bound  is used, a lower bound must also be used (e.g.,
      -INF; see below).  The lower and upper bounds  (or  if  just  the
      lower  bound  is  used,  then just this bound) may be omitted, in
      which case this form differs from the one described above only in
      so  far  as  with  this  form, the initial estimate and the FIXED
      attribute are enclosed in parentheses.  When FIXED is used, and a
      bound  appears  the bound must equal the initial estimate.  FIXED
      is implied when all three values are equal.  The lower bound  can
      be  -INF (i.e., -infinity), and the upper bound can be INF (i.e.,
      +infinity).  These are the defaults for lower and  upper  bounds.
      They  are  communicated  to NONMEM as numeric values -1000000 and
      1000000.

 ([low,] init [,up]) [FIXED]
      This is just like the form described above  except  that  if  the
      FIXED  attribute is used, the attribute occurs outside the paren-
      theses, and then if moreover, a bound appears, the bound need not
      equal the initial estimate.

 (low,,up)
      The commas are required.  Because no initial estimate is given, a
      search for an initial estimate (the Initial  Estimates  Step)  is
      undertaken  by  NONMEM.  With NONMEM 7.4, when initial thetas are
      to be estimated,  evaluations  can  now  be  done  for  FOCE  and
      Laplace, not just for FO.

 (value)xn
      Any  initial  value or group of initial values may be enclosed in |
      parentheses and followed by "xn", which means  to  replicate  the |
      values within parentheses n times ("repeated value").  The values |
      within the parenthesis may have any  of  the  above  forms.   For |
      example, the following two are equivalent:                        |

      $THETA 2 2 2 2 (0.001,0.1,1000) (0.001,0.1,1000) (0.001,0.1,1000) |
             (0.5 FIXED) (0.5 FIXED)                                    |

      $THETA (2)x4 (0.001,0.1,1000)x3 (0.5 FIXED)x2                     |

 UNINT (NM75)
      UNINT  is used during the Optimal Design Step to identify a theta
      as uninteresting.  UNINT may be used anywhere that FIXED  may  be
      used.

 label=value [FIXED] (NM75)
      This  is a compact method of defining an element of theta, speci-
      fying its initial estimate, and specifying a  symbolic  subscript
      for  this  element  of THETA.  The symbolic subscript may be used
      for THETA in abbreviated code, and will also identify  this  ele-
      ment  of THETA in the NONMEM output. (Only the first 9 characters
      of the label will appear).  If  new  $THETA  records  change  the
      ordering, the abbreviated code does  not have to be changed.  For
      example, suppose the third element of THETA that is defined  hap-
      pens to be
      $THETA CL=(0.0,7.0)
      The  abbreviated  code can use this symbolic subscript instead of
      the numeric subscript, e.g.,
      TVCL=THETA(CL)
      The NONMEM report will describe the relationship, e.g.,
      LABELS FOR THETAS THETA(3)=THETA(CL)
      and THETA(CL) rather than TH 3 will appear in the NONMEM report.

      It is also possible to specify the label using the
      $ABBR REPLACE THETA(CL)=THETA(3)
      control record but the NONMEM report will not identify the  rela-
      tionship  and  will  not generate the symbolic label THETA(CL) in
      the report.  For dynamic, or implicit mapping of labels, such  as
      for  various occasions, these still need to be done via the $ABBR
      REPLACE record.

      $THETA records must be placed ahead of any records that  use  the
      symbolic  label.  Therefore  these  records must be placed before
      $THETAI, $THETAR, and abbreviated code.

 NAMES (label ...)value ... (NM75)
      A compact way of defining one or more thetas with labels and ini-
      tial values. For example
      $THETA NAMES(V1,CL,Q,V2) (0.0,7.0) (0.0,7.0) (0.0,7.0) 7

 NUMBERPOINTS=n
      During  NONMEM's  search  for  an  initial  estimate, a number of
      points will be  examined.   This  number  will  be  automatically
      determined  by  NONMEM,  or it can be specified with this option.
      May also be coded NUMBERPTS, NUMPOINTS, NUMPTS.

 ABORT
      During the Initial Estimates step, NONMEM is to abort  when  PRED
      sets  the  error return code to 1.  (The PRED error return code n
      is set by the statement "EXIT n [k]" in abbreviated code,  or  by
      the  statement  IERPRD=n in user-supplied code, or by PREDPP when
      it detects an error.)  This is the default.

 NOABORT
      During the Initial Estimates step, NONMEM  is  simply  to  ignore
      values  of the theta vector that result in PRED error return code
      1.  (Ordinarily the first value of  the  theta  vector  is  never
      ignored.)  These will not be feasible values for an initial esti-
      mate.

 NOABORTFIRST
      Same as NOABORT option, but also applies to the  first  value  of
      the  theta  vector  that  is  tried.  It cannot be shortened; all
      characters must be coded.  May be used with the  NOABORT  option,
      in  which case the stronger condition (NOABORTFIRST) takes prece-
      dence.

 REFERENCES: Guide IV, section III.B.9 
