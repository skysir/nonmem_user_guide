


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                PRED                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PRED subroutine
 CONTEXT: User-supplied routine; required with NONMEM

 USAGE:

 Versions before NONMEM 7.2:

 SUBROUTINE PRED (ICALL,NEWIND,THETA,DATREC,INDXS,F,G,H)
 USE SIZES,     ONLY: DPSIZE,ISIZE,LVR
 REAL(KIND=DPSIZE) :: DATREC
 INTEGER(KIND=ISIZE) :: ICALL,NEWIND,INDXS
 DIMENSION :: THETA(*),DATREC(*),INDXS(*),G(LVR,*),H(LVR,*)

 With NONMEM 7.2 and higher:

 SUBROUTINE PRED (ICALL,NEWIND,THETA,DATREC,INDXS,F,G,H)
 USE SIZES, ONLY: DPSIZE,ISIZE
 USE PRDIMS,ONLY: GPRD,HPRD,GERD,HERD,GPKD
 IMPLICIT REAL(KIND=DPSIZE) (A-Z)
 REAL(KIND=DPSIZE) :: DATREC
 INTEGER(KIND=ISIZE) :: ICALL,NEWIND,INDXS
 REAL(KIND=DPSIZE) :: G(GPRD,*),H(HPRD,*)
 DIMENSION :: THETA(*),DATREC(*),INDXS(*)

 DISCUSSION:
 The  PRED  subroutine  is  called  by NONMEM to obtain modeled values.
 (PREDPP is a PRED subroutine that is distributed with NONMEM.)

 Input argument:

  ICALL

 ICALL=-1: the routine has been called for the  PRED_IGNORE_DATA   fea-
 ture  of  NONMEM  7.5.   One call per data record, at the start of the
 run.  These calls  occur  only  if  abbreviated  code  uses  variables
 PRED_IGNORE_DATA or PRED_IGNORE_DATA _TEST, or if the PRED_IGNORE_DATA
 option of $DATA is used.  Otherwise, there are no calls to  PRED  with
 ICALL=-1.

 ICALL=0:  PRED  has  been  called  for  initialization purposes at the
 beginning of the NONMEM run; one such call per run.   DATREC  contains
 the  first  data  record.  THETA contains the initial estimates.  PRED
 need not compute F, G, or H.

 ICALL=1: PRED has been  called  for  initialization  purposes  at  the
 beginning  of a NONMEM problem; one such call per problem.  Otherwise,
 identical to ICALL=0.

 ICALL=2: For the data record contained in DATREC, PRED has been called
 for  the  purpose  of computing F, the value of the prediction, and/or
 the values of other PRED-defined items, appropriate  for  the  record.
 PRED  should  compute  F, and also G and H as appropriate.  THETA con-
 tains the values to be used to compute F.   With  conditional  estima-
 tion, to obtain ETA values, PRED should call GETETA.

 ICALL=3:  PRED has been called for finalization purposes at the end of
 a NONMEM problem; one such call per (sub)problem.  DATEC contains  the
 first  data  record.   THETA contains the final estimates.  Otherwise,
 identical to ICALL=0.

 ICALL=4: For the data record contained in DATREC, PRED has been called
 during the Simulation Step for the purpose of computing a value of the
 dependent variable and, possibly,  values  of  independent  variables,
 appropriate  for  the record.  PRED should compute F (the value of the
 dependent variable).  THETA contains the initial estimates, which  are
 the  values  to  be  used  to compute F.  PRED should call SIMETA (and
 SIMEPS) to obtain simulated values of ETA (and EPS).

 ICALL=5: For the data record contained in DATREC, PRED has been called
 for the purpose of computing the expectation of the PRED item and pos-
 sibly, the expectations of other PRED-defined items,  appropriate  for
 the  record.  Such a call occurs when the marginal data item (MRG_) is
 defined in the data set and has a non-zero value for the  data  record
 in  question.   If the MRG_ data item on the record has the value 1 or
 2, the value returned by PRED in F contributes to the  expectation  of
 the  PRED  item.  Similarly, the values returned in other PRED-defined
 items contribute to expectations of these items.  THETA  contains  the
 final  estimates.  The expectations in question are over possible val-
 ues of ETA, and to obtain ETA values, PRED should call GETETA.

 ICALL=6:  PRED has been called for the purpose of  computing  the  raw
 data average of the DV data items and, possibly, the raw data averages
 of PRED-defined items.  Such a call occurs when the raw-data data item
 (RAW_)  is defined in the data set and has a non-zero value for a tem-
 plate data record (See template).  The value of the DV  data  item  in
 the  data  record contained in DATREC will be included in the raw data
 average of the DV data items.  However, when the raw data average cor-
 responding  to the label DV in a table or scatterplot is to be differ-
 ent from the raw data average of the DV  items  themselves,  PRED  may
 recompute  the value of DV.  PRED may return a value of 1 in F to omit
 the DV item in the data record from the average.

  NEWIND
      NEWIND=0: First record of the data set.  THETA value  may  differ
      from value at last call with this record.
      NEWIND=1: First record of the data set, THETA value does not dif-
      fer from value at last call with this record, and PRED is  nonre-
      cursive (see I_REC), or,
      First record of a subsequent individual record.
      NEWIND=2: Subsequent data record of an individual record.

  THETA
      The NONMEM THETA vector.

  DATREC
      The current data record.

  INDXS
      The  values specified in the $INDEX record of the NM-TRAN control
      stream.

 Output argument:

  F   When ICALL=2, the prediction associated  with  the  data  record.
      With  odd-type  data,  the  likelihood  of the observation in the
      record, but if there is no observation, F is ignored.

      When ICALL=4, the value of the simulated  observation  associated
      with  the  data record.  Alternatively, F can be ignored, and the
      DV item in the data record can be directly set to  the  value  of
      the  simulated  observation.   With  odd-type data, F is ignored;
      PRED should directly set the DV item to the value  of  the  simu-
      lated observation.

  G   An array of derivatives of F with respect to etas.  Values should
      be set when ICALL=2.
      G(i,1) is the partial derivative of F with respect to eta(i).
      When the data are population, G(i,j+1) is the second partial  de-
      rivative  of  F with respect to eta(i), eta(j) (lower-triangular;
      j=1,..., i).
      Second derivatives are needed only when the Laplacian  method  is
      used to estimate parameters.

  H   An  array  of  partial derivatives of F with respect to epsilons.
      When the data are population, values should be set when ICALL=2.
      H(i,1) is the derivative of F with respect to eps(i).
      H(i,j+1) is the partial derivative  of  H(i,1)  with  respect  to
      eta(j)
      These  mixed second derivatives are needed only when the INTERAC-
      TION option  is used to estimate parameters.  (See $ESTIMATION).

 Also see variables in NONMEM modules, NONMEM-PRED modules, and  PREDPP
 modules.
 (See variables in modules)

 REFERENCES: Guide I, section C.2 
 REFERENCES: Guide IV, section III.B.8 , IV 
 REFERENCES: Guide V, section 12.3 
