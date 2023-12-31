


 +--------------------------------------------------------------------+
 |                                                                    |
 |                 INITIALIZATION-FINALIZATION BLOCK                  |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for initialization and finalization
 CONTEXT: $PRED, $PK, $ERROR, $INFN abbreviated code

 SAMPLE:
 $PRED
     IF (ICALL.EQ.1) CALL SUPP(0,1)

 DISCUSSION:

 A  run  consists of one or more problems, and each problem consists of
 one or more subproblems, which are iterations of various tasks  speci-
 fied  in  the  problem. The "end" of a subproblem refers to the end of
 such an iteration.  Problems themselves may be organized  into  super-
 problems,  which  may  be  iterated.   The  "beginning" and "end" of a
 superproblem refers to the beginning and end of  the  first  and  last
 iterations  of the superproblem.  There are opportunities to make some
 rudimentary computations at:

 the beginning of a run (run initialization)
 the beginning of a superproblem (superproblem initialization)
 the beginning of a problem (problem initialization)
 the end of a subproblem (subproblem finalization)
 the end of a problem (problem finalization)
 the end of a superproblem (superproblem finalization)
 the end of a run (run finalization)

 For example, data transgeneration may take place at  problem  initial-
 ization.   Or, a variable may be initialized at problem initialization
 and modified at each subproblem  finalization,  and  its  final  value
 written  to a user file at problem finalization.  There is no opportu-
 nity to do subproblem initialization.  When  using  abbreviated  code,
 initialization  and finalization opportunities are signalled by values
 of the variable ICALL:

 ICALL=0
      Run initialization.

 ICALL=1
      Superproblem and Problem initialization.
      ICALL is set to 1 to signal problem  initialization.   When  this
      happens:

      If a superproblem requires initialization, test
      (i)    S1IT  (number  of current superproblem iteration) equals 1
      (or S2IT=1),
      (ii)  S1NUM (number of current superproblem)  equals  appropriate
      value, and
      (iii)  IPROB  (number  of current problem) equals number of first
      problem in superproblem,
      and if (i)-(iii) are true, do superproblem initialization
      (See Problem_Iteration_Counters).

 ICALL=3
      Subproblem, Problem, Superproblem and Run finalization.
      ICALL is set to 3 to signal problem finalization.  When this hap-
      pens:

      If a problem with subproblems requires finalization, test
      IREP=NREP  (number  of the current subproblem equals total number
      of subproblems) and if true, do problem finalization
      (See Simulation: NREP,IREP).

      If a superproblem requires finalization, test
      (i)  S1IT=S1NIT (number of  the  current  superproblem  iteration
      equals total number of superproblem iterations) (or S2IT=S2NIT),
      (ii)   S1NUM  (number of current superproblem) equals appropriate
      value,
      (iii) IPROB (number of current problem)  equals  number  of  last
      problem in the superproblem, and
      (iv)   IREP=NREP,  if there are subproblems with the last problem
      in the superproblem
      and if (i)-(iv) are true, do superproblem finalization.

      If a run having multiple problems requires finalization, test
      (i)  S1IT=S1NIT (or S2IT=S2NIT), if there are superproblems
      (ii)  S1NUM equals number of  last  superproblem,  if  there  are
      superproblems,
      (iii)  IPROB=NPROB  (number  of  the current problem equals total
      number of problems), and
      (iv)  IREP=NREP, if there are subproblems with the  last  problem
      in the superproblem
      and if (i)-(iv) are true, do run finalization.

 An  initialization  block  is a block of abbreviated code that is only
 executed at ICALL=0 or ICALL=1.  A finalization block is  a  block  of
 abbreviated code that is only executed at ICALL=3.  E.g.,

 IF (ICALL.EQ.1) THEN
  ... initialization block ...
 ENDIF

 Such  blocks may be present in $PRED, $PK, $ERROR, and $INFN blocks of
 abbreviated code.  If such blocks are present in  $PK  or  $ERROR,  an
 INFN routine is used to implement the logic in these blocks.  Initial-
 ization and finalization blocks will be implemented by means of a gen-
 erated FORTRAN subroutine.

 Variables  may  be  used as right-hand quantities even before they are
 defined; if an expression which  uses  such  a  variable  is  computed
 before  any  value of the variable is computed, the computation of the
 expression will be uncertain.

 Assignment, conditional, WRITE and PRINT statements may be used.

 In addition, these rules apply:

 Defined quantities are not regarded as random variables;  eta  deriva-
 tives are not computed in an initialization or finalization block

 Transgeneration  of  the  data is permitted.  NONMEM data items ID and |
 MDV may not be changed.  If a data item label may appear on  the  left
 of  an  assignment statement, then NM-TRAN generates assignment state-
 ments changing first the data item in the event or  data  record,  and
 then  the  value  of the local variable having that label.  Note, how-
 ever, that at ICALL=0,1,or 3, by default, references to data items are
 references  to those data items in the first data or event record.  To
 transgenerate an item in any  data  or  event  record  (including  the
 first),  use  of  the  NONMEM  utility  routine PASS is required.  See
 below.

 Calls to certain NONMEM routines are permitted:

      CALL PASS(MODE)
      CALL RANDOM(n,R)
      CALL SUPP(ie,ic)

      CALL PASS(MODE) must be coded  exactly  in  this  way.   If  CALL
      PASS(MODE)  is  present, MODE becomes a reserved variable and may
      be used only with other instances of CALL  PASS(MODE).   Multiple
      calls to PASS may be present.

      If CALL RANDOM(n,R) is present, R becomes a reserved variable and
      may be used only with other instances of CALL RANDOM(n,R).  n may
      only be an integer value 1-10.

      SUPP  is  used  to suppress portions of the NONMEM output report.
      (See supp).  The arguments ie and ic may only be 0 or 1.

 The following variables may be used on the right (their values  change
 with calls to PASS):

      Data record items (See PASS: PASSRC).
      NEWIND (See PASS NEWIND: NWIND).
      NIREC, NDREC (See Record Counters: NIREC,NDREC).
      NEWL2 (See PASS New L2 record: NEWL2).
      ETA when ICALL=3 (See Simulation: ETA,EPS)
      LIREC (See Size of Individual Record)
      PRED_, RES_, WRES_ when ICALL=3 (See PRED,RES,WRES).
      IERE, IERC when ICALL=3 (See Estim Covar Error Codes).
      NINDR, INDR1, INDR2 (See NINDR INDR1 INDR2).
      NREP, IREP (See Simulation: NREP,IREP).
      NPROB,  IPROB,  S1NUM, S2NUM, S1NIT, S2NIT, S1IT, S2IT (See Prob-
      lem_Iteration_Counters).

 RETURN and EXIT statements may be used.

 DOWHILE loops are permitted.  The syntax is as follows.

      DO WHILE (condition)
       .. statements ..
      END DO

      Here is an example of a transgeneration loop.

      $PRED
      IF (ICALL.EQ.0) THEN
           MODE=0
           CALL PASS (MODE)
           MODE=2
           CALL PASS (MODE)
           DO WHILE (MODE.EQ.2)
           ... transgeneration statements ...
           CALL PASS (MODE)
           ENDDO
          RETURN
      ENDIF

      This type of usage of the PASS routine can be coded more  simply,
      as follows:

      $PRED
      IF (ICALL.EQ.0) THEN
           DOWHILE (DATA)
           ... transgeneration statements ...
           ENDDO
          RETURN
      ENDIF

      The  DOWHILE  (DATA)  causes the transgeneration statements to be
      executed with each data record.  In effect, NM-TRAN supplies  the
      statements  MODE=...  and  CALL  PASS(MODE) that are shown in the
      above.

 Variables that are first defined in an initialization block or  final-
 ization  block  are  not  stored  globally  in NMPRD4, but rather, are
 stored in module PRINFN.  (This makes  them  available  to  subroutine
 MIX.)

 Variables  defined  in  an initialization or finalization block may be
 used freely outside of such blocks, and vice versa: PRED-defined vari-
 ables defined outside such blocks may be used within them.

 THETA variables may be used in initialization and finalization blocks.
 They are obtained from the subroutine argument.
 At ICALL=0, THETA contains the initial estimates for the  first  prob-
 lem.
 At ICALL=1, THETA contains the initial estimates for the current prob-
 lem.
 At ICALL=3, THETA contains the final estimates for the  current  prob-
 lem.

      Here is an example of code that could be used during a simulation
      with multiple subproblems.
      $PRED
      IF (ICALL.EQ.1) THEN
         SUM=0
         N=0
         RETURN
      ENDIF
      IF (ICALL.EQ.3) THEN
        N=N+1
        SUM=SUM+THETA(1)
        RETURN
      ENDIF
      IF (ICALL.EQ.3.AND.N.EQ.NREP) THEN
        MEAN=SUM/N
        PRINT *,MEAN
      ENDIF

 The following variables may be referenced as right-hand quantities  in
 initialization and finalization blocks.  (Note also that unsubscripted
 arrays may appear in a WRITE statement.)

      OMEGA(n,m)
      SIGMA(n,m)
      (See Parameter Values: Initial and Final).
      (See write print).

 The following variables may be referenced as right-hand quantities  in
 finalization blocks.
 (Note also that unsubscripted arrays may appear in a WRITE statement.)

      SETHET(n)
      SETHETR(n)
      SEOMEG(n,m)
      SESIGM(n,m)
      OBJECT
      (See Standard Errors).
      (See Objective Function Value).
      (See write print).

 The following statements are forbidden:

      CALL SIMETA(ETA)
      CALL SIMEPS(EPS)

 When the PRED repetition feature is used, the variables RPTO and PRDFL
 may appear as left-hand quantities in initialization blocks.
 (See Repetition Variables).

 (See abbreviated).

 REFERENCES: Guide II, section D.2.2 
 REFERENCES: Guide VI, section VI.A , Figure 37
