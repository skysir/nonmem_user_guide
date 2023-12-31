


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           COPYING BLOCK                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code especially for copying pass
 CONTEXT: $PRED, $PK, $ERROR, $AES, $DES abbreviated code

 SAMPLE:
  IF (COMACT.EQ.1) TVCL=CL

 DISCUSSION:

 Values  of variables displayed in tables and scatterplots are obtained
 (i.e. copied) from module NMPRD4.  There  are  particular  times  when
 data  records  are  passed  to PRED for the purpose of obtaining these
 values; these are called copying passes.  The variable COMACT  signals
 that  a copying pass is in progress when its value is positive.  There
 may be a number of copying passes, but if  (and  only  if)  values  of
 variables  are  to  be  displayed, there is at least one copying pass.
 With the first copying pass, the value of COMACT  is  1  and  that  of
 MIXNUM  is  1.  If a mixture model with k subpopulations is used, then
 COMACT remains 1 during a total of k copying  passes,  and  with  each
 copying  pass  the  value  of  MIXNUM increments by 1.  If conditional
 estimation is used (or if the POSTHOC option appears), then there fol-
 lows a set of copying passes where the value of COMACT is 2 and again,
 MIXNUM increases from 1 to k.  If a mixture model is not  used,  there
 are  altogether  at  most  two  copying  passes, one with COMACT=1 and
 another with COMACT=2.

 If values of a variable from earlier  copying  passes  are  needed  in
 later passes, the values for the variable should be stored in the SAVE
 region of module NMPRD4.  (See comsav).  When the values are stored in
 the SAVE region, that value computed with a given data record during a
 copying pass will be found in NMPRD4 when the same  record  is  passed
 during  the  next  copying pass, i.e. it will have been saved from the
 previous copying pass.  This is in contrast  to  the  usual  behaviour
 (with a noncopying pass), where with a given data record, the value in
 NMPRD4 is the value computed with the previous data record.  The value
 is 0 when a record is passed during the first copying pass, and though
 PRED may set it or reset it during a  subsequent  copying  pass,  this
 need not be done (see discussion below about COMACT).

 The  values used in tables and scatterplots (whether or not these val-
 ues are stored in the SAVE region) are those copied from  NMPRD4  with
 the last copying pass.

 A  copying  block is a block of abbreviated code that is only executed
 during a copying pass, i.e. when COMACT has a positive value.  Special
 rules apply, which allow the user to be less concerned about using the
 COMSAV option of the $ABBREVIATED record (see example below):

 COM(i) variables that are defined in a copying block are  referred  to
 as  explicit  SAVE variables.  PRED-defined variables that are defined
 in a copying block, other than COM(i) variables, are  referred  to  as
 implicit  SAVE  variables.   Collectively, the two types of SAVE vari-
 ables are referred to as SAVE variables.  NM-TRAN  stores  SAVE  vari-
 ables in the SAVE region of module NMPRD4.

 Implicit and explicit SAVE variables cannot both appear in abbreviated
 code.  The COMRES option of the $ABBREV record cannot be used when any
 implicit  SAVE  variables  are used, but it must be used when explicit
 SAVE variables appear (as whenever COM(i) variables appear).

 The size of the SAVE region of NMPRD4 depends on the COMSAV option  of
 the $ABBREV record.  This option may be used in three ways:

      No COMSAV option.  The SAVE region of module NMPRD4 will nonethe-
      less include all the SAVE variables.

      COMSAV=n (n>=0).  NM-TRAN will, if necessary, extend the size  of
      the  SAVE  region from n to a larger value so that all SAVE vari-
      ables will be included in the region.

      COMSAV=-1.  There is to be no SAVE region.  Variables defined  in
      a copying block will not be SAVE variables.

 EXAMPLE OF USAGE:

 $ERROR
 Y=F+F*EPS(1)
 IPRED=F
 IF (COMACT.EQ.1) FT=F
 WR=(DV-IPRED)/FT

 $EST ... POSTHOC
 $TABLE FT IPRED WR

 With  the  first  copying pass the value of COMACT is 1, which signals
 that during this copying pass, all ETA variables are set to 0.   Since
 the  option  POSTHOC appears, with a subsequent copying pass the value
 of COMACT is 2, which signals that during this copying pass,  all  ETA
 variables are set to their conditional estimates.

 In  this example, WR is set to the weighted intra-individual residual.
 When COMACT=1, the prediction for the typical  individual  (ETA=0)  is
 computed  and stored in FT.  FT is a SAVE variable, so this value will
 have been saved, and when COMACT=2, this same value will be  found  in
 FT.   With  COMACT=2, the weight used with the residual is computed to
 be the reciprocal of this FT value, while IPRED is computed from  con-
 ditional  estimates of the ETA variables and thus its value applies to
 to the individual with these estimates, rather  than  to  the  typical
 individual.   The  values  of  IPRED and WR appearing in the table are
 those obtained during the last copying pass.    The  tabled  value  of
 IPRED is based on the conditional estimate of ETA.  The value of FT is
 also the value obtained during the last copying pass, but it  in  turn
 is  also  the value obtained from the first copying pass, as the value
 of FT has not changed during any subsequent  copying  pass,  and  this
 value is based on ETA=0.

 Since  FT  is  a SAVE variable, a SAVE region will have been allocated
 where this variable will be stored, and (unless there is  some  reason
 to  use  the  COMSAV option other than to insure this) the user not be
 concerned about this option.

 (See COMACT,COMSAV)
 (See PRED-Defined Variables).
 (See abbreviated).

 REFERENCES: None.
