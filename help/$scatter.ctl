


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            $SCATTERPLOT                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Requests that NONMEM generate one or more scatterplots
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SCATTERPLOT  list1 VS list2 [BY list3]
               [FROM n1] [TO n2] [UNIT]
               [ORD0|NOORD0] [ABS0|NOABS0] [FIRSTONLY] [OBSONLY]
               [NOSUB=[0|1]]
               [UNCONDITIONAL|CONDITIONAL] [OMITTED]

 SAMPLE:
 $SCATTERPLOT      (RES WRES) VS TIME BY ID

 DISCUSSION:
 Requests  that  families of NONMEM scatterplots be produced.  Up to 20
 families of scatterplots may be included in the problem.  May also  be
 coded $SCATTERS or $SCATTERGRAMS.

 OPTIONS:

 list1
      A  list  of  item  labels to be plotted on the ordinate axis (the
      long axis on the printed output).  The list may  be  enclosed  in
      parentheses. It may include:

      Data item labels.

      Special items PRED, RES, and WRES.

 The  user  may  request  the  following additional diagnostic items by
 including their name in the list.

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

 Labels ETA(1), ETA(2), ... ,ETA(10), ... , ETA(70), etc.,
 or alternatively, labels ETA1, ETA2, ... ,ET10, ... , ET70, etc., cor-
 responding to eta(1), eta(2), etc.
 The labels in the output will be ETA1, ETA2, ...  ,  ET10,  ...  ET70,
 etc.

      With NONMEM 7.3, a range of etas may be requested:                |
      $SCAT ETAS(1:2) VS ETA3                                           |
      is equivalent to                                                  |
      $SCAT ETA1 ETA2 VS ETA3                                           |
      However,  unlike  $TABLE, $SCAT will ignore implied endings, such |
      as                                                                |
      $SCAT ETAS(1:LAST) VS ETA3                                        |
      And just interpret it as                                          |
      $SCAT ETA1 VS ETA3                                                |
      With NONMEM 7.4, more flexible syntax is available, using TO  and |
      BY.                                                               |
      (See $table).                                                     |

      With  NONMEM 7.4, a symbolic label specified in $ABBR REPLACE may |
      be listed in $SCAT. For example:                                  |
      $ABBR REPLACE ETA(CL)=ETA(1)                                      |
       ...                                                              |
      $SCAT  ETA(CL) VS ETA3                                            |

 Reserved positions of MODULE  NMPRD4  (see  $ABBREV).   COM(k)  or  :k
 denotes  the  kth  reserved position.  (There must be exactly 4 digits
 after ":". Use leading 0's as necessary.)  E.g., COM(3) or :003.
 Labels of the form :k will be used in the output.

 Labels of PRED-defined items in MODULE NMPRD4 if abbreviated  code  is
 present  (up  to  PDT distinct such labels in any one problem, for all
 tables and scatterplots.  PDT is  a  constant  in  resource/SIZES.f90;
 default  value  is  500.)   These  may  include labels of the NM-TRAN-
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

 Synonyms may be defined on either the $TABLE or  $SCATTER  record  for
 special  items  PRED,  RES, WRES; for PRED-defined items; for NM-TRAN-
 defined items; and for reserved positions of MODULE NMPRD4.

 E.g., assume that IWRES is a PRED-defined label, that at least 3 posi-
 tions of NMPRD4 are reserved, and that NM-TRAN has generated A00032 as
 the label for a derivative in the  generated  FSUBS  routine.   Either
 $TABLE or $SCATTER records may include:

         WRES=RES1,IWRES=RES2,COM(3)=ABC,0032=DK.

 For  a  discussion of the values of ETAs and PRED-defined items (e.g.,
 are they based on initial or final values of theta?  Simulated or zero
 or conditional values of eta?), see values.

 Elements of G and H
 E.g., $SCATTER G11 BY G21
 The  format  is  Gk1  or Hk1, where k is an integer value, e.g. 1-9 or
 01-99 or 001-999.  Gk1 requests the value of G(k,1), and Hk1  requests
 the  value  of H(k,1), where G and H are arguments of subroutine PRED.
 G(k,1) is the partial of F (the prediction) with  respect  to  ETA(1),
 and  H(k,1)  is  the  partial of F with respect to EPS(1).  (HH may be
 coded instead of H, but it is treated as if it were H.)  A variable of
 the  form  Gk1  or Hk1 is not a reserved variable. If it is previously
 defined (i.e., if it is listed in $INPUT,  or  used  on  the  left  in
 abbreviated  code, or used as a synonym e.g., $TABLE G11=COM(1)), then
 that definition of the variable is used, and there is no  change  from
 previous versions of NM-TRAN.  Only if there is no other previous def-
 inition of the variable will it be understood to be an element of G or
 H.
 What  NM-TRAN  actually displays is the variable in MODULE NMPRD4 that
 was used to compute the derivative of  interest  (a  variable  in  the
 series  A00nnn, C00nn, or D00nnn) with the appropriate synonym such as
 G11.  If there is no such variable, this is an error. NMTRAN will  not
 display  variables  that  are  not computed , e.g., G41 when there are
 only 3 etas in the problem, or when there are 4 etas but  ETA(4)  does
 not contribute to the value of Y.
 [There  is a workaround if the zero is wanted as a place holder in the
 table.  In abbreviated code ($ERROR or $PRED or $PK)
 G41=0 .

      (See Displayed PRED-Defined Items).

 list2
      Like list1, but includes the labels of items to be plotted on the
      abscissa  axis  (the short axis on the printed output).  The list
      may be enclosed in parentheses.  The word "VS"  is  optional;  it
      may  be  omitted  if the lists are enclosed in parenthesis, or if
      each list consists of exactly one item  label.  VS  may  also  be
      coded *.  Each pair of labels, one from list1 and one from list2,
      defines a family of scatterplots.

 list3
      A list of one or two item  labels  or  synonyms.   Each  pair  of
      labels  from  list1  and list2 produces a family of scatterplots,
      one scatterplot for each unique value (or combination of  values)
      of  the data item(s) in list3.  If the BY option is omitted, each
      pair of labels from list1 and list2 produces a  family  comprised
      of a single scatterplot.

 UNIT A line of unit slope is superimposed on the scatterplots.

 ORD0 A  line  through  the zero value on the ordinate axis is superim-
      posed on the scatterplots.  May also be coded ORDZERO.

 NOORD0
      Prevents a zero line from being superimposed on the ordinate axis
      of the scatterplots.  May also be coded NOORDZERO.
      If  neither  ORD0  nor  NOORD0  is  present, NONMEM automatically
      superimposes a zero line on the ordinate axis if it is  appropri-
      ate for the type of data item.

 ABS0 A  zero line is superimposed on the abscissa axis of the scatter-
      plots.  May also be coded ABSZERO, AB0, or ABZERO.

 NOABS0
      Prevents a zero line from being superimposed on the abscissa axis
      of  the  scatterplots.   May  also  be  coded  NOABSZERO, NOABS0,
      NOABZERO, NOAB0.
      If neither of ABS0 and NOABS0 is  present,  NONMEM  automatically
      superimposes  a zero line on the abscissa axis if it is appropri-
      ate for the type of data item.

 FIRSTONLY
      Only the first data record from each individual record  may  con-
      tribute   a   point  to  the  scatterplot.   May  also  be  coded
      FIRSTRECORDONLY or FIRSTRECONLY.

 OBSONLY
      The scatterplot will only use  data  records  with  MDV=0.   This
      option  applies  independently of FIRSTONLY.  It is not necessary
      when either DV, RES, or WRES is plotted.

 FROM n1
      n1 is the number of the first data record which may  "contribute"
      to the scatterplot.  Default: n1 is 1.

 TO n2
      n2  is  the number of the last data record which may "contribute"
      to the scatterplot.  Default:  There  is  no  upper  limit.   All |
      appropriate  records  will  contribute.  To restore the NONMEM VI |
      behavior, use TO  n1+899.                                         |

 With FIRSTONLY, n1 and n2  refer  to  (first  records  of)  individual
 records.  The remaining options apply to all $SCATTERPLOT records.

 NOSUB=[0|1]
      With  NOSUB=0,  label substitution will be performed scatterpots.
      (See $ABBREVIATED).  This is the default.   With  NOSUB=1,  label
      substitution will not be performed.

 CONDITIONAL
      The Scatterplot Step is implemented only when the Estimation Step
      terminates successfully or  is  not  implemented.   This  is  the
      default.

 UNCONDITIONAL
      The Scatterplot Step is always implemented.  This is the default.

 OMITTED
      The Scatterplot Step is not implemented.

 When  DV,  RES, or WRES is plotted, records having MDV=1 are not plot-
 ted.

 The following symbols are used in scatterplots:

               *   1 point
   Overstriking:
             2-9   2-9 points
             A-Z   10-35 points (A=10, B=11, ... , Z=35)
               $   more than 35 points

 REFERENCES: Guide IV, section III.B.17 
 REFERENCES: Guide V, section 9.5.2 , 10.7.2 
