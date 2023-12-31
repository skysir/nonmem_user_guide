


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         DATA_AVERAGE BLOCK                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for computation of raw-data-averages
 CONTEXT: Abbreviated code

 SAMPLE:
  $ERROR
 IF (ICALL.EQ.6) THEN
  ... data average block ...
 ENDIF

 DISCUSSION:
 A  data average block is a block of abbreviated code that is only exe-
 cuted when ICALL=6.  This value of ICALL  occurs  when  the  raw-data-
 average (RAW_) data item is defined in the data set and has a non-zero
 value for some records.  Data average blocks are not required when the
 raw-data-average  data  item is present, but they allow the user addi-
 tional functionality.  Such blocks may be present in $PRED and $ERROR.
 If  for  a given observation record matching a template record, a data
 average block sets the DV variable to a value different from  the  one
 in the record, this value is the one included in the average.

 In  this  example,  the displayed DV value for each template record is
 the proportion of DV items greater than 10 on those  records  matching
 the template record.

  Y=...
  TRDV=DV
  IF (ICALL.EQ.6) THEN
     DV=0
     IF (TRDV.GT.10) DV=1
  ENDIF

 PRED  may return a value of 1 in F, indicating that the DV item in the
 record is not to be included in the  average.   Continuing  the  above
 example,  suppose  that the DV variable is set to a value only when it
 exceeds 2, so that the  average  is  the  proportion  of  observations
 exceeding  10 among those that exceed 2.  The code is as follows.  See
 note 2 below on the returned value of F.

 Y=...
 TRDV=DV
 IF (ICALL.EQ.6) THEN
    IF (TRDV.LE.2) Y=1
    DV=0
    IF (TRDV.GT.10) DV=1
 ENDIF

 When ICALL=6, PRED and ERROR are called  with  successive  records  as
 usual.   However,  with each observation record, all output from these
 routines other than a value set for the DV data item  and  values  set
 for  PRED-defined  (ERROR-defined) items V used in a table or scatter-
 plot and located in the SAVE region (See save_region (and the F value,
 in  so  far as it is or is not 1) is ignored.  In the same way that an
 average is formed for the DV, averages are formed for the elements  of
 V.   Upon  entry  into PRED or ERROR with a given record, the value of
 the DV item is the one on the record, and the values  for  V  are  the
 values  that  will  be used in a table and/or scatterplot.  If PRED or
 ERROR changes one of these values, the new value is the  one  used  in
 the  average.   If a value is not changed, the unchanged value is used
 in the average.  (For a nonobservation record, the  output  from  that
 record is completely ignored.)

 The  following  series  of  examples  concern taking averages of PRED-
 defined items.  They involve a mixture model, where, under each of the
 subpopulations, there is a parameter PA whose value depends on an eta.
 There is a template record.

 Example A.

 Suppose first that PA does not depend  on  (interindividually-varying)
 covariate  values.  Suppose moreover, that during copying passes, with
 each individual record a value of the quantity  Q=P1*PA1+P2*PA2+P3*PA3
 has  been  computed  and stored in the SAVE region.  Here, the P's are
 the mixture probabilities, and the PA's are conditional  estimates  of
 PA  under  the different models for the different subpopulations.  Due
 to the presence of the template record, at ICALL=6, the average  of  Q
 across  the  eta  estimates  (from individuals with the same covariate
 values as are contained in the template record)  will  be  computed  -
 without any need for a data average block.  More precisely, at ICALL=6
 a pass through the data set occurs, during  which  a  value  of  Q  is
 obtained  with each of the observation records that match the template
 record.  However, the average of the Q values is an average of within-
 individual  averages,  and  since Q does not vary within an individual
 record, the within-individual average computed for that record is  the
 same value of Q as is obtained with each of the observation records of
 the individual record (matching the template record).   The  resulting
 average  of  the  Q  values is an estimate of the expected value of PA
 over the subpopulations and the randomly-varying PA's.

 Example B.

 Suppose that PA depends on a covariate X
 E.g. PA=THETA(1)*X**THETA(2)*EXP(ETA(1)),
 and that one is interested in an estimate of the expected value of  PA
 for  X=x.  Suppose also that the value for X in the template record is
 x.  Once again, suppose that during copying passes, with each individ-
 ual  record  a  value  of  Q  has been computed and stored in the SAVE
 region - using whatever value of X appears in the  individual  record.
 At  ICALL=6, if and only if the observation records within an individ-
 ual record have the value X=x  and  match  the  template  record  with
 respect  to  the other relevant data items, will the within-individual
 average be included in the average  Q.   Moreover,  these  observation
 records are the very ones whose value for Q is of interest.

 Example C.

 However,  under a well-specifed model, ETA(1) should be independent of
 X, and so one might want to use an average of Q across the  eta  esti-
 mates from all individuals (with observation records).  Then one might
 (i) during copying passes, for each individual record compute Q  using
 the  specific value X=x - regardless of what value of X appears in the
 individual record, and (ii) include the record

 $OMIT X

 to prevent the values of X from affecting the match with the  template
 record.

 Example D.

 The  strategy  in example C fails if it is necessary to match on X for
 the purpose of forming averages other than the average Q, and it is at
 least  awkward  if one is interested in the expected value of PA for a
 variety of values of X.  Here is an alternative strategy.

 $ABB COMRES=4 COMSAV=4
 $ERROR
  ...
 Y=...
 IF (COMACT.EQ.2) THEN
    IF (MIXNUM.EQ.1) COM(1)=ETA(1)
    IF (MIXNUM.EQ.2) COM(2)=ETA(2)
    IF (MIXNUM.EQ.3) COM(3)=ETA(3)
 ENDIF
 IF (ICALL.EQ.6) THEN
    PA1=THETA(1)*TEMPLT(X)**THETA(2)*EXP(COM(1))
    PA2=THETA(1)*TEMPLT(X)**THETA(2)*EXP(COM(2))
    PA3=THETA(1)*TEMPLT(X)**THETA(2)*EXP(COM(3))
    COM(4)=MIXP(1)*PA1+MIXP(2)*PA2+MIXP(3)*PA3
 ENDIF

 $TABLE COM(1) COM(2) COM(3) NOPRINT FILE=junk
 $TABLE ID ... COM(4) ...

 During the copying passes, the values of the eta estimates are  stored
 (in  items  COM(1),  COM(2) and COM(3) of the SAVE region) rather than
 the value of Q.  The first table record  appears  because  unless  the
 items  COM(1), COM(2), and COM(3) are displayed, at ICALL=6 their val-
 ues are 0.  TEMPLT(X) refers generically to the value of X on the tem-
 plate record (see Special Rule 4 below); so there can be numerous tem-
 plate records with different values of X.

 Example E.

 Examples A-C are not explicit about how the mixture probabilities  are
 computed.   They  might be obtained via MIXP, in which case during the
 copying passes (similar to what happens at ICALL=6  with  example  D),
 with  a given individual record, they are the probabilities pertaining
 to that record.  As long as the mixture probabilities do not depend on
 interindividual-varying  covariates, the mixture probabilities are the
 same no matter what individual record it is  to  which  they  pertain.
 But if the probabilities depend on interindividual-varying covariates,
 and especially if one wants to estimate the expected value  of  Q  for
 numerous different sets of values for these covariates, then one might
 use:

 $ERROR
 include nonmem_reserved_general                                        |
 Y=...
 IF (COMACT.EQ.2) THEN
    IF (MIXNUM.EQ.1) COM(1)=ETA(1)
    IF (MIXNUM.EQ.2) COM(2)=ETA(2)
    IF (MIXNUM.EQ.3) COM(3)=ETA(3)
 ENDIF
 IF (ICALL.EQ.6) THEN
    PA1=THETA(1)*TEMPLT(X)**THETA(2)*EXP(COM(1))
    PA2=THETA(1)*TEMPLT(X)**THETA(2)*EXP(COM(2))
    PA3=THETA(1)*TEMPLT(X)**THETA(2)*EXP(COM(3))
    COM(4)=MIXPT(1)*PA1+MIXPT(2)*PA2+MIXPT(3)*PA3                       |
 ENDIF

 The mixture probabilities found in MIXPT  pertain  to  the  individual
 record  containing  the  template  record  (see Special Rule 7 below).
 Because they have been computed using the  covariate  values  in  that
 record,  they  are  commensurate with the way PA has been computed for
 the different subpopulations.  With versions prior to NONMEM 7.3,  the |
 probabilities  in MIXPT must be referenced via verbatim code (note the |
 double quotes).                                                        |

 $ERROR                                                                 |
 " USE ROCM_REAL MIXPT=>MIXP_RAW                                        |
  ...                                                                   |
 "  COM(4)=MIXPT(1)*PA1+MIXPT(2)*PA2+MIXPT(3)*PA3                       |

 Special rules apply to data average blocks and passes through the data
 with ICALL=6.

 1)   No eta derivatives are computed in a data average block.

 2)   A  RETURN  statement  may be used in a data average block.  If Y=
      appears (i.e. Y is assigned a value)  in  a  data  average  block
      before  the  RETURN (not necessarily in the same block containing
      the RETURN), then F is set to Y (F=Y); otherwise F is set to to 0
      (F=0).   If there is no RETURN statement in a data average block,
      then as usual, F is set to the value of Y assigned  by  the  time
      PRED (ERROR) exits.

 3)   Loops  are  permitted  in a data average block.  The syntax is as
      follows.

      DO WHILE (condition)
       .. statements ..
      END DO

 4)   During calls with ICALL=6, the template data record is  found  in
      NONMEM global variable TEMPLT.  Items in the template data record
      may be referred to in abbreviated code of the data average  block
      by position or by label, e.g., TEMPLT(1) or TEMPLT(ID).

 5)   During  calls  with  ICALL=6,  the  repetition feature may not be
      used.

 6)   If a mixture model is used, then during calls with ICALL=6,  both
      MIXNUM  and  MIXEST are the index of the subpopulation into which
      the individual (whose data record is being passed) has been clas-
      sified.

 7)   If a mixture model is being used, then during calls with ICALL=6,
      (the final estimates of)  the  mixture  probabilities  associated
      with  the  individual  record  containing the template record are
      found in NONMEM global variable  MIXPT.   Code  that  uses  these
      probabilities must be verbatim code.

 (See raw_).
 (See MIX CONTR: TEMPLT)
 (See Mixture model: MIXPT).

 REFERENCES: None.
