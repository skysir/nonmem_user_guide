


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          SIMULATION BLOCK                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for simulation
 CONTEXT: Abbreviated code

 SAMPLE:
 $PK
  IF (ICALL.EQ.4) CL=THETA(1)+ETA(1)

 DISCUSSION:
 A  "simulation block" is a block of abbreviated code that is only exe-
 cuted when ICALL=4 (during simulation).  Such blocks may be present in
 $PK,  $ERROR,  and $PRED, and may be implemented by means of generated
 FORTRAN subroutines.  E.g.,

 IF (ICALL.EQ.4) THEN
  ... simulation block ...
 ENDIF

 Special rules apply to such blocks.

 1)   No eta derivatives are computed in a simulation block.

 2)   Transgeneration is permitted.  NM-TRAN allows a data  item  label
      to appear on the left of an assignment statement.  NM-TRAN gener-
      ates assignment statements changing first the data  item  in  the
      event  or  data  record,  and then the local variable having that
      label.  E.g.,  suppose WT is listed in $INPUT:

      IF (ICALL.EQ.4) WT=70+70*ETA(3)

      The generated code is:
              IF(ICALL.EQ.4.D0)THEN
              EVTREC(NVNT,8 )=70.D0+70.D0*ETA(03)
              WT=EVTREC(NVNT,08)
              ENDIF

      NONMEM and PREDPP reserved data items should not be modified dur- |
      ing  simulation.   Transgeneration  is  permitted with simulation |
      with subproblems.  With all versions of NONMEM, the data set  for |
      each  subproblem after the first is the same data set used by the |
      previous subproblem, and includes any  changes  (transgeneration) |
      made by the previous subproblem.

 3)   Calls to certain NONMEM routines are permitted:

      CALL SIMETA(ETA)
      CALL SIMEPS(EPS)
      CALL RANDOM(n,R)
      where   n  is  an  integer  1-10.   If  CALL RANDOM is present, R
      becomes a reserved variable used for the random number.

      Note that NM-TRAN provides the  necessary  calls  to  SIMETA  and
      SIMEPS  in generated routines.  Explicit calls are used in abbre-
      viated code only to obtain different values of ETA and EPS.

 4)   A RETURN statement may be used.  If in $ERROR or $PRED,  and  the
      RETURN  occurs  in  a  simulation block, then Y may be assigned a
      value prior to the return.  If so, then F is set (F=Y); otherwise
      F is not set.

 5)   Loops are permitted.  The syntax is as follows.

      DO WHILE (condition)
       .. statements ..
      END DO

      Here are some examples.  For a truncated normal distribution that
      only requires testing the eta value directly, this  code  can  be
      used:
      $PK
       ...
      DO WHILE (ETA(1).GT.5)
      CALL SIMETA(ETA)
      ENDDO

      Another example:
      IF (ICALL.EQ.4.AND.NEWIND.NE.2) THEN
        DO WHILE (ETA(1).GT..5.OR.ETA(1).LT.-.5)
          CALL SIMETA(ETA)
        ENDDO
      ENDIF
      IF (ICALL.EQ.4) WT=70+70*ETA(1)

      (With  these  two examples, the first random seed of the $SIMULA-
      TION record must have the NEW option.  Note also that, because of
      the  previous  automatic  call to SIMETA, ETA(1) requires no ini-
      tialization, but that R in the next example does.)

      IF (ICALL.EQ.4.AND.NEWIND.NE.2) THEN
        R=1
        DO WHILE (R.GT..5.OR.R.LT.-.5)
          CALL RANDOM(2,R)
        ENDDO
      ENDIF
      IF (ICALL.EQ.4) WT=70+70*R

      This example illustrates how a categorical variable  with  equal-
      likely  probabilities  can  be  generated from a random number R,
      uniformly distributed between 0 and 1.  In this example, the cat-
      egorical variable BIN takes values 1 through 5.

      IF (ICALL.EQ.4) THEN
        CALL RANDOM(2,R)
        BIN=INT(R*5)+1
      ENDIF

      The number 5 can be replaced with any other positive integer n to
      obtain an n-valued categorical variable.  Here INT is  the  func-
      tion  that  transforms  a  nonnegative number x into the greatest
      integer not exceeding x.  The effect of this simulation  code  is
      to perform the transformation:

      BIN=1 if R < .2
      BIN=2 if R < .4 and R >= .2
      BIN=3 if R < .6 and R >= .4
      BIN=4 if R < .8 and R >= .6
      BIN=5 if R <  1 and R >= .8

 6)   The "EXIT n k" statement may be used.  The value of n may be 0, 1
      or 2.  The value of k is referred to as the PRED EXIT CODE.
      If it is desired that the simulation be  immediately  terminated,
      then use an EXIT 2 code:
      IF(ICALL==4.and.IPRED<0.1 .and. TIME>20.0) EXIT 2

      With  versions  of NONMEM prior to 7.2, the "EXIT 1" statement in
      the Simulation step also caused NONMEM to abort.   As  of  NONMEM
      7.2, if an error occurs in PREDPP during simulation such as
      PK PARAMETER FOR KA IS NON-POSITIVE
      or  a  user-implemented  EXIT 1 is issued during simulation, then
      PRED will be called with a new ETA  and  EPS.   This  feature  is
      referred  to  as  Simulation Error Forgiveness.  NONMEM describes
      this as PRED SIMULATION REDO  in  the  NONMEM  report  file.   It
      writes to the NONMEM report file a description of the data record
      and THETA and ETA values, for example
       PRED SIMULATION REDO PRED EXIT CODE = 1 INDIVIDUAL NO.       1
       ID= 1.00000000000000E+00   (WITHIN-INDIVIDUAL) DATA REC NO.   1
       THETA=
        3.00E+00   8.00E-02   4.00E-02
       ETA=
        4.66E-01   2.91E-03   9.95E-01
       MESSAGE ISSUED FROM SIMULATION STEP

      If ten such errors occur in the same subject, then it is supposed
      that  the  cause  of  the simulation error is not due to an occa-
      sional bad random sample, but is caused by a systematic error  in
      the  control  stream file. The simulation step is terminated with
      the message
      PRED ERROR OCCURRED TOO OFTEN ON SIMULATION
      instead of a message
      SIMULATION STEP PERFORMED

      With NONMEM 7.5, the PRED  EXIT  CODE  k  may  be  in  the  range
      1000-9999. For example,
      IF(ICALL==4.and.IPRED<0.01 .and. TIME>20.0) EXIT 1 2300
      This can only occur with user's EXIT code; PREDPP will not gener-
      ate this kind of EXIT.  NONMEM will try PRED SIMULATION  REDO  up
      to  10000  times.   The  message "PRED SIMULATION REDO" itself is
      written to PRDERR  up to 30 times.   After  that,  the  following
      message is written to PRDERR:
       SUBSEQUENT PRED SIMULATION REDO ERROR MESSSAGES SUPPRESSED
      NONMEM continues trying new ETA and EPS. Be careful that the con-
      dition does not occur too often (causing  wasteful  computation).
      After  10000  tries, the simulation is terminated as a protection
      against an infinite loop. The following  message  is  written  to
      PRDERR:
       TOO MANY CONSECUTIVE PRED ERRORS (>10000) OCCURRED ON SIMULATION

 (See PRED Exit Code).
 (See abbreviated).
 (See INTRODUCTION TO NONMEM 7, Simulation Error Forgiveness (NM72)).
 (See INTRODUCTION TO NONMEM 7, Extensions to Simulation Error Forgive-
 ness (NM75)).

 REFERENCES: Guide IV, section III.B.13 , IV.I 
 REFERENCES: Guide V, section 12.4.8 
 REFERENCES: Guide VI, section III.E.2 , IV.B.2 
