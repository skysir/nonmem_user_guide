


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                INFN                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: INFN subroutine
 CONTEXT: User-supplied subroutine; replaces a PREDPP dummy routine.

 USAGE:
 SUBROUTINE INFN (ICALL,THETA,DATREC,INDXS,NEWIND)
 USE SIZES, ONLY: ISIZE,DPSIZE
 USE NMPRD_INT, ONLY: NWIND
 INTEGER(KIND=ISIZE) :: ICALL,INDXS(*),NEWIND
 REAL(KIND=DPSIZE)   :: THETA(*),DATREC(*)

 DISCUSSION:
 A  run  consists of one or more problems, and each problem consists of
 one or more subproblems.  With NONMEM 75, there is an  opportunity  to
 drop  data  records  from the NONMEM data set at the start of the run,
 using the PRED_IGNORE_DATA feature.  There are opportunities  to  make
 some  rudimentary computations at the beginning of a run (run initial-
 ization), at the beginning of a problem (problem  initialization),  at
 the  end  of  a  subproblem (subproblem finalization), at the end of a
 problem (problem finalization), and at the end of a run (run finaliza-
 tion).   For  example, at problem initialization, data transgeneration
 may take place, or a variable that  will  be  modified  at  subproblem
 finalization,  may be initialized.  At problem finalization, the value
 of this variable may be written to a user file.  There is no  opportu-
 nity to do subproblem initialization.

 The  INFN  subroutine  is  called by PREDPP to make initialization and
 finalization computations.  The version distributed with PREDPP  is  a
 "stub"  or  "dummy" routine that does nothing. It may be replaced by a
 user-written code.

 Input/Output argument:

  ICALL

      ICALL=-1: INFN may  set  PRED_IGNORE_DATA=1  for  records  to  be
      dropped.  The following is required:
           USE NMPRD_INT, ONLY: PRED_IGNORE_DATA,PRED_IGNORE_DATA_TEST
      The  option $DATA ... PRED_IGNORE_DATA  is required to cause NON-
      MEM to make a PRED_IGNORE_DATA  pass.

      ICALL=0: INFN may now make computations for  run  initialization.
      ICALL may be reset by INFN to a number in the range 1-8999.  This
      number will appear on NONMEM output, allowing the user  to  iden-
      tify the INFN routine being used.

      ICALL=1:  INFN  may now make computations for problem initializa-
      tion.

      ICALL=3: INFN may now make computations for subproblem  finaliza-
      tion.

      ICALL=3: INFN may now make computations for problem finalization.
      If there are subproblems, first do  subproblem  finalization,  if
      required.   Then  test  for  IREP=NREP (the number of the current
      subproblem equals the total number of subproblems), and if  true,
      do problem finalization.  Values of IREP and NREP may be found in
      NONMEM modules
      (See Simulation: NREP,IREP).

      ICALL=3: INFN may now make  computations  for  run  finalization.
      First  do  problem  finalization,  if  required.   Then  test for
      IPROB=NPROB (the number of the current problem equals  the  total
      number of problems), and if true, do run finalization.  Values of
      IPROB and NPROB may be found in NONMEM modules
      (See Problem_Iteration_Counters).

 Input arguments:

  DATREC
      DATREC contains the first data record  of  the  current  problem.
      (At  ICALL=0,  the  current problem is the first problem.)  Using
      PASS (see below), the contents of DATREC are  replaced  by  other
      data  records  for  the current problem, allowing all these other
      records to be read and even modified.

  THETA
      The NONMEM THETA vector.
      At ICALL=0, THETA contains the initial estimates  for  the  first
      problem.
      At  ICALL=1, THETA contains the initial estimates for the current
      problem.
      At ICALL=3, THETA contains the final estimates  for  the  current
      problem.

  INDXS
      The  values specified in the $INDEX record of the NM-TRAN control
      stream.  (This is the NONMEM INDXS array starting at position 12,
      the first position beyond those positions used by PREDPP itself.)

  NWIND
      NWIND has value 0 when INFN is called.  It changes value during a
      pass through the data using PASS (see below).
      NWIND=0: First record of the data set.
      NWIND=1: First record of a subsequent individual record.
      NWIND=2: Subsequent data record of an individual record.

 EXAMPLES OF USAGE:

 Initialization
      Constants in the PK routine can be set if they are  stored  in  a
      module or common block declared both in INFN and PK.

 Transgeneration
      The  data can be accessed and even modified via use of the NONMEM
      utility routine PASS in routine INFN  (See pass).   As  the  data
      records  are  passed one-by-one to INFN, each record is stored in
      turn in DATREC.
      Data can be transgenerated and additional data items can be  pro-
      duced  at  both initialization and finalization.  Data items that
      appear in a table or scatterplot with a given subproblem or prob-
      lem  are  produced or left unchanged at the subproblem or problem
      finalization.
      At finalization estimates of the eta's can be  obtained  via  the
      NONMEM  utility  routine  GETETA.   When used in conjunction with
      PASS, the values returned for the eta's with each call to  GETETA
      are appropriate for the individual whose data record is currently
      in DATREC.

 Interpolation
      As an  example  of  transgeneration,  interpolated  values  of  a
      covariate  can  be computed for event records in which values are
      missing,  e.g.  for  other-type  event  records  that  have  been
      included  so  that  predictions  can  be obtained at the times in
      these records. This could also be done in PK or ERROR,  but  then
      this  would be done with every call to these routines; if done in
      INFN, the computation is done once only.

 (See Infn interpolation example).

 REFERENCES: Guide VI, section VI.A , Figure 37
