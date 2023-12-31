


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $TABLE                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Requests that NONMEM generate a table
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $TABLE  [list1] [BY list2]
         [PRINT|NOPRINT] [FILE=filename]
         [NOHEADER|ONEHEADER] [ONEHEADERALL]
         [NOTITLE|NOLABEL]
         [FIRSTONLY|LASTONLY|FIRSTLASTONLY] [NOFORWARD|FORWARD]
         [APPEND|NOAPPEND]
         [FORMAT=s] [LFORMAT=s] [RFORMAT=s]
         [IDFORMAT=s]
         [NOSUB=[0|1]]
         [EXCLUDE_BY list3]
         [PARAFILE=[filename|ON|OFF]
         [ESAMPLE=n1][WRESCHOL]
         [SEED=n2] [CLOCKSEED=[0|1]]
         [RANMETHOD=[n|S|m] ]
         [VARCALC=[0|1|2|3]]
         [FIXEDETAS=(list)]
         [NPDTYPE=[0|1]]
         [UNCONDITIONAL|CONDITIONAL] [OMITTED]

 SAMPLE:
 $TABLE          ID DOSE WT TIME

 DISCUSSION:
 Requests that a NONMEM table be produced.  Up to 10 $TABLE records may
 be included in a given problem.

 OPTIONS:

 list1
      A list of item labels (i.e., user-chosen item types) to be tabled
      along with DV and the special items PRED, RES, and WRES.

           The  user may request the following additional special diag-
           nostic items by including their name in the list.

      NPRED, NRES, NWRES
           Calculated assuming non-conditional estimation and  no  eta-
           epsilon  interaction.   NPRED  and NRES are same as PRED and
           RES. NWRES is same as WRES when INTERACTION is  not  set  in
           $EST.

      PREDI, RESI, WRESI
           Calculated assuming  non-conditional  estimation  with  eta-
           epsilon interaction.  Always same as PRED, RES, and WRES.

      CPRED, CRES, CWRES
           Calculated  assuming  conditional  estimation  and  no  eta-
           epsilon interaction.

      CPREDI, CRESI, CWRESI
           Calculated  assuming conditional estimation with eta-epsilon
           interaction.

      CIPRED, CIRES,CIWRES
           Conditional individual values

      CIPREDI, CIRESI,CIWRESI
           Conditional individual values with eta-epsilon  interaction.

      EPRED, ERES, EWRES
           Monte-Carlo  generated  diagnostics  and  are not linearized
           approximations like the other diagnostic types. EWRES is the
           Monte-Carlo version of CWRESI.

      ECWRES
           Monte-Carlo version of CWRES.

      NPDE
           Monte-Carlo generated  normalized  probability  distribution
           error.

      NPD  The correlated (or non-decorrelated) NPDE value.

      OBJI
           Objective function values for each individual (same as given
           in the root.phi file).  The sum of the individual  objective
           function values is equal to the total objective function.

 In addition the list may include:

 Data item labels.

 Labels ETA(1), ETA(2), ... ,ETA(10), ... , ETA(70), etc.,
 or  alternatively,  labels  ETA1, ETA2, ... ,ETA10, ... , ETA70, etc.,
 corresponding to eta(1), eta(2), etc.
 The labels in the output will be ETA1, ETA2, ...  ,  ET10,  ...  ET70,
 etc.
      With NONMEM 7.3, a range of etas may be requested:                |
      ETAS(k:n)                                                         |
      is equivalent to                                                  |
      ETAk, ..., ETAn                                                   |
      where  n  >  k.  LAST can be used in place of n, and requests the |
      last (highest numbered) eta in the problem.  E.g. ETAS(1:LAST)    |
      With NONMEM 7.4, more flexible syntax is available:               |
      The word TO may be used in place of ":".                          |
      The BY expression may be used:                                    |
      ETAS(1 TO 10 by 3) prints out etas 1,4,7,10                       |
      ETAS(LAST TO 1 by -3) prints out etas 10,7,4,1 (assuming LAST=10) |
      A number list may be given:                                       |
      ETAS(1,5,12,4) prints out etas 1, 5, 12, 4.                       |
      ETAS(4:1) prints etas 4, 3, 2, 1                                  |
      ETAS(4:1 by -2) prints etas 4, 2                                  |
      ETAS(1:4 by -1) prints etas 4, 3, 2, 1 (the by value sets the direction).|

      With NONMEM 7.4, a symbolic label specified in $ABBR REPLACE  may |
      be listed in $TABLE. For example:                                 |
      $ABBR REPLACE ETA(CL)=ETA(1)                                      |
       ...                                                              |
      $TABLE ETA(CL)                                                    |

 Reserved  positions of MODULE NMPRD4 (See $ABBREVIATED).  COM(k) or :k
 denotes the kth reserved position.  (There must be  exactly  4  digits
 after ":". Use leading 0's as necessary.)  E.g., COM(3) or :003.
 Labels of the form :k will be used in the output.

 Labels  of  PRED-defined items in MODULE NMPRD4 if abbreviated code is
 present (up to PDT distinct such labels in any one  problem,  for  all
 tables  and  scatterplots.   PDT  is a constant in resource/SIZES.f90;
 default value is 500.)  These  may  include  labels  of  the  NM-TRAN-
 defined items:
   0nnn     e.g., 0010   stands for A00nnn
   1nnn     e.g., 1010   stands for A01nnn
   2nnn     e.g., 2010   stands for C00nnn
   3nnn     e.g., 3010   stands for D00nnn
   4nnn     e.g., 4010   stands for E00nnn
   5nnn     e.g., 5010   stands for F00nnn
   6nnn     e.g., 6010   stands for P00nnn

 These may also include:
   labels  VECTRA(1),  VECTRA(2),  ...  ,VECTRA(9),  or  alternatively,
   labels VA_1, VA_2, ... ,VA_9, corresponding to VECTRA(1), VECTRA(2),
   ..., VECTRA(9).
   The labels in the output will be VA_1, VA_2, ... , VA_9.
   Similarly, for VECTRB and VECTRC.

 Synonyms  may  be  defined on either the $TABLE or $SCATTER record for
 special items PRED, RES, WRES; for PRED-defined  items;  for  NM-TRAN-
 defined items; and for reserved positions of MODULE NMPRD4.

 E.g., assume that IWRES is a PRED-defined label, that at least 3 posi-
 tions of NMPRD4 are reserved, and that NM-TRAN has generated A00032 as
 the  label  for  a  derivative in the generated FSUBS routine.  Either
 $TABLE or $SCATTER records may include:

         WRES=RES1,IWRES=RES2,COM(3)=ABC,0032=DK.

 For a discussion of the values of ETAs and PRED-defined  items  (e.g.,
 are they based on initial or final values of theta?  Simulated or zero
 or conditional values of eta?), see values.

 Elements of G and H
 E.g., $TABLE G11 G21 G31 H11 H21
 The format is Gk1 or Hk1, where k is an integer  value,  e.g.  1-9  or
 01-99  or 001-999.  Gk1 requests the value of G(k,1), and Hk1 requests
 the value of H(k,1), where G and H are arguments of  subroutine  PRED.
 G(k,1)  is  the  partial of F (the prediction) with respect to ETA(1),
 and H(k,1) is the partial of F with respect to  EPS(1).   (HH  may  be
 coded instead of H, but it is treated as if it were H.)  A variable of
 the form Gk1 or Hk1 is not a reserved variable. If  it  is  previously
 defined  (i.e.,  if  it  is  listed  in $INPUT, or used on the left in
 abbreviated code, or used as a synonym e.g., $TABLE G11=COM(1)),  then
 that  definition  of the variable is used, and there is no change from
 previous versions of NM-TRAN.  Only if there is no other previous def-
 inition of the variable will it be understood to be an element of G or
 H.
 What NM-TRAN actually displays is the variable in MODULE  NMPRD4  that
 was  used  to  compute  the  derivative of interest (a variable in the
 series A00nnn, C00nn, or D00nnn) with the appropriate synonym such  as
 G11.   If there is no such variable, this is an error. NMTRAN will not
 display variables that are not computed , e.g.,  G41  when  there  are
 only  3  etas in the problem, or when there are 4 etas but ETA(4) does
 not contribute to the value of Y.
 [There is a workaround if the zero is wanted as a place holder in  the
 table.  In abbreviated code ($ERROR or $PRED or $PK)
 G41=0.
 Now G41 may be listed in $TABLE or $SCATTER.]
 This  feature  is  designed  so  that  the  verbatim code in the "com-
 pute.cwres" R documentation is unnecessary.
 E.g., Instead of:
   $ABB COMRES=5
   "LAST
   "  COM(1)=G(1,1)
   "  COM(2)=G(2,1)
   "  COM(3)=G(3,1)
   "  COM(4)=HH(1,1)  (or H(1,1) with $PRED)
   "  COM(5)=HH(2,1)  (or H(2,1) with $PRED)
 $TABLE ID COM(1)=G11 COM(2)=G21 COM(3)=G31
  COM(4)=H11 COM(4)=H21
  IPRED MDV NOPRINT ONEHEADER FILE=cwtab1
 Use only:
 $TABLE ID G11 G21 G31 H11 H21
  IPRED MDV NOPRINT ONEHEADER FILE=cwtab1

      (See Displayed PRED-Defined Items).

      When tables are printed, the maximum number of  labels  permitted
      in  list1  is  8;  otherwise,  it  is PDT.  (But see the NOAPPEND
      option.)

 list2
      A list comprised of one or more labels from list1.  The  rows  of
      the  table  are sorted on the data items in list2.  List2 may not
      appear when the number of labels in  list1  is  greater  than  8.
      That  is,  a  table  with  more than 8 data items also may not be
      sorted.

 NOSUB=[0|1]
      With NOSUB=0, label substitution  will  be  performed  for  final
      estimates  in  table  files.   (See $ABBREVIATED).   This  is the
      default.  With NOSUB=1, label substitution will not be performed.

 list3
      A list comprised of one or  more  items  that  are  permitted  in
      list1, e.g., data item labels and labels of PRED-defined items in
      MODULE NMPRD4.  They follow option EXCLUDE_BY.  Labels  in  list3
      are  called  exclusion  variables.  If one or more of them have a
      non-zero value for a given data record, the row of the table cor-
      responding  to  the  data  record will be excluded from the table
      file.  Exclusion variables are not  listed  in  the  table  file.
      They  have no effect on the printed table or scatters in the NON-
      MEM output, e.g., they do not cause any rows to be  deleted  from
      the printed table and are displayed in the printed table.

 PARAFILE=filename
      Weighted  residuals are evaluated before the first $TABLE record.
      As of NONMEM 7.4, this  computation  is  parallelized  if  paral-
      lelization is on when the first Table Step is implemented.
      PARAFILE=filename  specifies  a  different parafile than was used
      for the previous step.
      PARAFILE=ON turns on parallelization for the weighted residuals.
      PARAFILE=OFF turns off parallelization for the  weighted  residu-
      als.
      The  PARAFILE  option  may be specified on any $TABLE record, but
      applies to all $TABLE records.

 ESAMPLE=n1
      n1 specifies the number of random samples used to  calculate  the
      Monte-Carlo  diagnostics.   Should be specified only on the first
      $TABLE record.  Default is 300.

 WRESCHOL (NM73)                                                        |
      Use the Cholesky square root of  the  variance,  rather than  the |
      eigenvalue   square  root,  when  computing  weighted  residuals. |
      Should be specified only on the first $TABLE  record.   This  can |
      speed up the Table Step when there are a large number of observa- |
      tions per individual.

 SEED=n2
      n2 specifies the starting seed  for  generating  the  Monte-Carlo
      diagnostics.   Should  be  specified  only  on  the  first $TABLE
      record.  Default is 11456.

 CLOCKSEED=[0|1] (NM75)
      If CLOCKSEED=1 (default is  0),  actual  starting  seed  will  be
      10000*(seconds  after  midnight)+SEED.   This  allows  a  control
      stream to produce  different  stochastic  results  for  automated
      replications,  without  the  need to modify the seed value in the
      control stream file in each replication.

 RANMETHOD=[n|S|m]
      n: the random number generator used for  the  Monte-Carlo  simua-
      tions of weighted residual items.
        0: ran0 of reference [5], minimal standard generator
        1: ran1 of reference [5], Bays and Durham.
        2: ran2 of reference [5].
        3: ran3 of reference [5], Knuth. (Default)
        4: NONMEM's traditional random number generator used in $SIMULATION
      S: sobol sequence                                                 |
      m: the type of scrambling desired                                 |
        0: no scrambing (S0 is the same as S)                           |
        1: Owen type scrambling                                         |
        2: Faure-Tezuka type scrambling                                 |
        3: Owen plus Faure-Tezuka type scrambling.                      |
      See  the  description  of  RANMETHOD for $ESTIM.  Among the Sobol |
      sequence methods, the S2 method  appears  to  provide  the  least |
      biased  random samples, that is nearly uniform distribution, with |
      good mixing in multi-dimensional spaces.
      See INTRODUCTION TO NONMEM 7, Reference [5]
      See INTRODUCTION TO NONMEM 7, Monte Carlo Importance Sampling EM
      RANMETHOD should be specified only on the first  $TABLE  command.
      The  RANMETHOD  set  in  the $TABLE command does not propagate to
      $EST or $CHAIN.

      Options PRINT, NOPRINT, HEADER, NOHEADER, NOLABEL, NOTITLE, FILE,
      FIRSTONLY,  LASTONLY,FIRSTLASTONLY,  FORWARD,  NOFORWARD, APPEND,
      NOAPPEND, FORMAT, VARCALC  apply to the individual $TABLE record.
      They must be specified for each table to which they apply.

 PRINT
      A  printed  table  appears  in  the  NONMEM  output.  This is the
      default.

 NOPRINT
      No printed table appears in the NONMEM output.

 FILE=filename
      The table is written to the given file in character  form,  e.g.,
      ASCII  or  EBCDIC,  according to the hardware platform.  Filename
      may not contain embedded spaces.  If it  contains  commas,  semi-
      colons,  or  parentheses,  then  it  must be surrounded by single
      quotes ' or double quotes ".  Filename  may  also  contain  equal
      signs  if it is enclosed in quotes.  Filename may contain at most
      71 characters.  If filename is the same as any option of the $TA-
      BLE  record,  it must be enclosed in quotes.  Filename can differ
      between $TABLE records.
      Default: No table file is output.  Required with NOPRINT.

 NOHEADER
      Used only with the FILE option.  No header lines are included  in
      the table file.

 ONEHEADER
      Used only with the FILE option.  Only the first line of the table
      is a header line.

 ONEHEADERALL (NM74)
      Used only with the FILE option and FORWARD.  Only the first  line
      of the table file is a header line.  May also be coded ONEHEADER-
      PERFILE.

 NOLABEL
      Used only with the FILE option.  Do not print column labels.   It |
      may  be  combined  with  ONEHEADER to print only the title at the |
      beginning of each table.                                          |

 NOTITLE                                                                |
      Used only with the FILE option.  Do not print table  titles.   It |
      may be combined with ONEHEADER to print only the column labels at |
      the beginning of each table.  NOLABEL NOTITLE  is  equivalent  to |
      NOHEADER.

 FIRSTONLY
      Only information corresponding to the first data record from each
      individual record appears  in  the  table.   May  also  be  coded
      FIRSTRECORDONLY or FIRSTRECONLY.

 LASTONLY
      Only  information corresponding to the last data record from each
      individual record appears in the table.  May also be  coded  LAS-
      TRECORDONLY or LASTRECONLY.

 FIRSTLASTONLY
      Only  information corresponding to the first and last data record
      from each individual record appears in the table.

 NOFORWARD
      Used only with the FILE option.  When the table  file  is  opened
      during a given (sub)problem, it is positioned at the start of the
      file.   This is the default.  However, when  there  are  multiple
      $TABLE  records within the same problem and having the same file-
      name, the situation is a little more complicated;  see  the  text
      describing the FORWARD option.

 FORWARD
      Used only with the FILE option.  When a table file is opened dur-
      ing a given (sub)problem, it is forwarded to the end of the file. |
      This  allows a table file to accumulate tables from multiple sub- |
      problems and superproblems.  Moreover, if in the same  (sub)prob-
      lem  the  $TABLE record is followed by a contiguous succession of
      additional $TABLE records having the same filename as  the  given
      record,  then  even  though  some of these additional records may
      have the NOFORWARD option, or have neither the NOFORWARD nor  the
      FORWARD  options,  the  FORWARD  option will apply to each of the
      records in the succession.

 APPEND
      Items DV, PRED, RES, WRES appear automatically as the last 4 col-
      umns of the table.  This is the default.

 NOAPPEND
      Requests that items DV, PRED, RES, WRES not appear automatically.
      When this is specified, the number of labels  (i.e.,  user-chosen
      item  types)  that  may appear in the table can be as large as 12
      (rather than 8) for a  printed  table,  and  as  large  as  PDT=4
      (rather  than  PDT) for a table file.  If items PRED, RES, and/or
      WRES are explicitly coded in list1, then they appear in  the  ta-
      ble, exactly as listed.  (Previously to NONMEM VI 2.0, they could
      be included in the list, but were suppressed from the portion  of
      the table described by list1 in favor of the automatically-gener-
      ated items.)

 FORMAT=s
      This option defines the delimiter and  number  format  for  table
      files.  It affects table files until a different FORMAT is speci-
      fied.  s defines the delimiter [,|s(pace)|t(ab)]  followed  by  a
      Fortran  format  specification.   The default is s1PE11.4.  There
      are many more options for FORMAT.
      For more details, see the format help item:
      (See format).

      Alternately,  use  LFORMAT and/or RFORMAT.
      See INTRODUCTION TO NONMEM 7, FORMAT=s1PE11.4

 LFORMAT=s
      Specifies the format of the full label record of a table.  Allows
      different formats for different columns.  Sample:
      LFORMAT="(4X,A4,4(',',4X,A8))"

 RFORMAT=s
      Specifies  the  format  of  the  full  numeric record of a table.
      Allows different formats for different columns.  Sample:
      RFORMAT="(F8.0,,4(',',1PE12.5))"

      Multiple LFORMAT options and RFORMAT options may be  present  and
      will  be concatenated.  The format specifications are not checked
      by NMTRAN. If either is invalid, the run  will  fail  in  NONMEM.
      Both  LFORMAT  and  RFORMAT  affect table files until a different
      format is specified.  Use  LFORMAT="NONE"  or  RFORMAT="NONE"  to
      resume  use  of the default format (which may have been specified
      by the FORMAT option) in a subsequent problem.                    |

 IDFORMAT=s(NM75)                                                       |
      By default the ID column has the same format as specified by FOR- |
      MAT.   IDFORMAT  specifies  the format for the ID column in table |
      files.  If an improper format is given, it defaults  to  that  of |
      FORMAT.  Some examples:                                           |

      IDFORMAT=I                                                        |
           Integer value, left adjusted in the field.                   |

      IDFORMAT=I6                                                       |
           Integer  value,  right adjusted in the first 6 characters of |
           the field.                                                   |

      IDFORMAT=F6.1                                                     |
           Floating value, with single digit to the right of the  deci- |
           mal.                                                         |

 VARCALC=[0|1|2|3](NM74,NM75)                                           |
      To report standard errors associated with etas (individual param- |
      eters) in the  tables  for  user-defined  variables,  set  $TABLE |
      ...VARCALC=1.   See  setest.ctl  in the examples directory.  This |
      appends an item named item_SE following each user-defined item in |
      list1.   If  using  RFORMAT formatting option, make sure to allow |
      enough format fields to include  reported  standard  errors.   In |
      addition, full variances-covariances among all user-defined vari- |
      ables and PREDPP parameters will be outputted  to  file  root.vpd |
      (the FORMAT used for this file is that defined in the $EST state- |
      ment).  With NONMEM 7.5, file root.vpt  is  also  created.   This |
      contains  variance-covariances  associated with thetas (or omegas |
      and sigmas) as well  as  those  associated  with  etas/individual |
      parameters.   To  append  the  comparable  total  standard errors |
      item_SE to tables, set $TABLE ...VARCALC=3.  To only  create  the |
      vpd  and  vpt  files,  and not report SE's to the table, set VAR- |
      CALC=2.  This option must be re-coded for each $TABLE record  for |
      which  SE's  are  wanted.  If VARCALC=1 or 2 or 3 is requested at |
      least once among any of the tables, then the variance  items  are |
      written  to  the  vpd  and vpt files.  VARCALC=0 requests neither |
      SE's nor vpd nor vpt, and is  the  default.   Values  of  COMACT, |
      which  identify copying passes, may need to be tested in abbrevi- |
      ated code when this feature is used.                              |
      See INTRODUCTION TO NONMEM 7, Requesting Standard Errors to User- |
      Defined and PREDPP Parameters                                     |

 FIXEDETAS=(list)(NM74)                                                 |
      Specified  etas  may be treated as if they are fixed effects when |
      evaluating population diagnostics during the $TABLE  step.   This |
      is  particularly  suitable  for  super-ID  $LEVEL  etas that span |
      groups of subjects, as if they were a fixed effect when  evaluat- |
      ing  populations  characteristics during the $TABLE step, such as |
      PRED, CWRES, NPDE, etc. In this way, the PRED evaluated will  be, |
      not  of  the total population, but of a given site level for that |
      subject.  List is a number-list of etas.   For  example,  FIXEDE- |
      TAS=(3-6)                                                         |
      A  number-list  may contain a single integer, a range of integers |
      (with -), or a series of integers and ranges separated by comma.  |

 NPDTYPE=1                                                              |
      The strict stochastic (Monte Carlo) method over the data y domain |
      as well as etas is implemented for NPD diagnostics.               |

 NPDTYPE=0                                                              |
      An  asymptotic  assessment  of  the residual variability is used. |
      This is the default.

 UNCONDITIONAL
      The Table Step is always implemented.  This is the default.

 CONDITIONAL
      The Table Step is implemented only when the Estimation Step  ter-
      minates successfully or is not implemented.

 OMITTED
      The Table Step is not implemented.

 REFERENCES: Guide IV, section III.B.16 
 REFERENCES: Guide V, section 9.5.1 , 10.7.1 
