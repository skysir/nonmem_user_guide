


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               ERROR                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: ERROR subroutine
 CONTEXT: User-supplied subroutine; for use with PREDPP

 USAGE: Versions before NONMEM 7.2:

 SUBROUTINE ERROR (ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,F,G,HH)
 USE SIZES, ONLY: DPSIZE,ISIZE,LVR
 INTEGER (KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS
 REAL(KIND=DPSIZE) :: THETA,F,G,HH
 DIMENSION IDEF(*),THETA(*),EVTREC(IREV,*),INDXS(*),G(LVR,*)
 DIMENSION HH(LVR,*)

 With NONMEM 7.2 and higher:

 SUBROUTINE ERROR (ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,F,G,HH)
 USE SIZES,     ONLY: DPSIZE,ISIZE
 USE PRDIMS,    ONLY: GERD,HERD
 IMPLICIT REAL(KIND=DPSIZE) (A-Z)
 REAL(KIND=DPSIZE) :: EVTREC
 INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS
 DIMENSION :: IDEF(*),THETA(*),EVTREC(IREV,*),INDXS(*)
 REAL(KIND=DPSIZE) :: G(GERD,*),HH(HERD,*)

 DISCUSSION:
 The  ERROR  subroutine  is  called by PREDPP to model intra-individual
 error in observed values.  It can also be used to to  convert  predic-
 tions  from  PREDPP, i.e., scaled drug amounts, to other types of pre-
 dictions (for example, to obtain the prediction of a drug effect as  a
 function of concentration, in a pharmacodynamic study).

 Input argument:

  ICALL
      ICALL=1:  ERROR  has been called for initialization at the begin-
      ning of a NONMEM problem; one such call per problem.  EVTREC con-
      tains  the  first event record.  THETA contains the initial esti-
      mates.  ERROR must return values in IDEF which inform PREDPP what
      tasks it will perform at later calls.  It may also set elments of
      HH to the appropriate derivatives, when  it  has  requested  that
      ERROR be called only once per problem.

      ICALL=2:  ERROR has been called to obtain the modeled value; mul-
      tiple calls occur.  ERROR should compute derivatives of F for HH.
      If  ERROR  changes F, it should also change the derivates of F in
      G.

      ICALL=4: ERROR has been called during the Simulation Step; multi-
      ple calls occur.  ERROR should compute simulated observations and
      store them in F.

      ICALL=5: ERROR has been called during the computation of expecta-
      tions;  multiple calls occur.  Such calls occur when the marginal
      (MRG_) data item is defined in the data set and has non-zero val-
      ues  for  some records.  If the MRG_ data item has the value 1 or
      2, expectations are computed  for  that  record,  and  the  value
      returned  by  ERROR  in  F contributes to the expectation that is
      being computed for the PRED data item.  The value  of  an  ERROR-
      defined  item  contributes  to  the  expectation computed for the
      item.

      ICALL=6: ERROR has been called during the computation of raw data
      averages;  multiple  calls occur.  Such calls occur when the raw-
      data (RAW_) data item is defined in the  data  set  and  has  the
      value 1 for some records (See template).  ERROR may re-compute DV
      when the value of DV in the data record is not the quantity to be
      included  in  the average.  ERROR may return a value of 1 in F if
      no DV item is to be included in the average with  the  particular
      record.

  THETA
      The NONMEM THETA vector.

  EVTREC
      The PREDPP event record.

  INDXS
      The  values specified in the $INDEX record of the NM-TRAN control
      stream.  (This is the NONMEM INDXS array starting at position 12,
      the first position beyond those positions used by PREDPP itself.)

 Output argument:

  IDEF
      ERROR  should store values in IDEF only when ICALL=1.  A value of
      1 in IDEF(1) indicates that ERROR will store the  derivatives  of
      log  y  in HH (which causes PREDPP to multiply each element of HH
      by F).  That is, PREDPP understands that the "error in log y"  is
      modeled.

      The value in IDEF(2) describes when ERROR should be called:
      -1 call with every event record.
       0 call once per observation record.
       1 call once per individual record.
       2 call once per problem.

      The  value in IDEF(3) describes whether ERROR uses derivatives of
      compartment amounts (i.e. whether compartment amounts  themselves
      are used as random variables in arithmetic statements in ERROR).
      -1 ERROR may use derivatives of A.
       0 ERROR does not use derivatives A.
       1 ERROR does use derivatives of A.

      The  default  used  by PREDPP is IDEF(3)=-1.  However, when ERROR
      does not use A, then if IDEF(3) is set to  0,  PREDPP  can  avoid
      some  time-consuming processing.  Indeed, when $ERROR abbreviated
      code is supplied, and there is  no  reference  to  a  compartment
      amount  A(n) (as a random variable in an arithmetic statement) in
      the abbreviated code (or to its derivatives  in  verbatim  code),
      then NM-TRAN sets IDEF(3)=0.

 Input/Output argument:

  F   On  input,  the  prediction  based  on the pharmacokinetic model,
      i.e., the value of the scaled drug amount in the observation com-
      partment.   On output, F may be unchanged, or it may be modified,
      e.g., when a PD prediction is needed.  If ERROR  modifies  F  and
      uses  population  eta  variables  to do so, then at ICALL=2 ERROR
      must call GETETA to obtain eta values prior to modifying F.  When
      ICALL=4,  ERROR should calculate the simulated observation (after
      calling  SIMETA  and/or  SIMEPS  to  obtain  simulated  etas  and
      epsilons, as necessary) and place its value in F.  With the Simu-
      lation Step, ERROR may return the simulated observation as the DV
      data item, rather than in the argument F.  With odd-type data the
      simulated observation must be  returned  as  the  DV  data  item.
      (See $estimation).

  G   On  input,  at  ICALL=2 when the data are population, an array of
      partial derivatives of F with respect to etas
      G (i,1) is the partial derivative of F with respect to eta(i).
      G (i,j+1) is the second derivative of F with respect  to  eta(i),
      eta(j) (lower-triangular; j=1, ..., i).
      (Second derivatives are only needed with estimation by the Lapla-
      cian method.)
      On input, G contains zeros when the data are single-subject data.
      At ICALL=2, ERROR must modify G when it has  changed  F  and  has
      thus  changed the derivatives of F with respect to the population
      etas.

  HH  An array of partial derivatives of F with respect to  etas  (when
      the  data are single-subject data) or epsilons (when the data are
      population).  Values should be stored when ICALL=2 and also  when
      ICALL=1 if ERROR sets IDEF(2)=2.
      HH(i,1) is the derivative of F with respect to eta(i) or eps(i).
      HH(i,j+1) is the derivative of H(i,1) with respect to eta(j) (but
      are only needed with conditional estimation when  the  dependence
      on  etas  of the variance of intra-individual random error should
      be preserved in the computation of the  objective  function;  see
      the INTERACTION option (See $ESTIMATION).

 Also  see  NONMEM  read-only  modules (of the form ROCMn), NONMEM-PRED
 modules (of the form NMPRDn), and PREDPP  read-only  modules  (of  the
 form PROCMn).

 REFERENCES: Guide IV, section V.C.6 
 REFERENCES: Guide V, section 8 
 REFERENCES: Guide VI, section IV 
