


 +--------------------------------------------------------------------+
 |                                                                    |
 |                  COMPARTMENT INITIALIZATION BLOCK                  |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for compartment initialization
 CONTEXT: $PK abbreviated code

 SAMPLE:
 $PK
 IF (A_0FLG.EQ.1) THEN
  ... compartment initialization block ...
 ENDIF

 DISCUSSION:
 A  "compartment  initialization  block" is a block of abbreviated code
 that sets the initial state of the kinetic system.  It is to  be  exe-
 cuted  only  when  A_0FLG=1.  Such a block may be present only in $PK.
 PREDPP sets A_0FLG to 1 at a call to PK with the first event record of
 an individual record (if the data are population data), with the first
 event record of the data set (if the data  are  single-subject  data),
 and with a reset record.

 Special rules apply to such a block.

 1)   Values  may be assigned to reserved variables A_0(n), but only in
      a compartment initialization block.  The value of the  amount  in
      the  nth compartment (the nth element of the state vector) is set
      to the value assigned to A_0(n).   If PK is called  with  a  dose
      record  or  a  dose-reset record where the dose is input into the
      nth compartment, this amount is then increased by the  amount  of
      the  (bioavailable) dose.  If a value is assigned to A_0(n), then
      it is not necessary that values be assigned to any of the remain-
      ing  variables  A_0(m).  A value to the output compartment cannot
      be assigned.  A_INITIAL(n) is a synonym for A_0(n).

 2)   The statement "IF (A_0FLG.EQ.1)" and  the  corresponding  "ENDIF"
      statement  may  be  included explicitly in abbreviated code, thus
      defining an  explicit  compartment  initialization  block.   (See
      example 1 below.)

 3)   A_0(n)  may  be assigned a value with an unconditional statement.
      This defines an implicit compartment initialization block; NMTRAN
      inserts  "IF  (A_0FLG.EQ.1) ..." before the statement and "ENDIF"
      after it.  (See example 2 below.)
      Indicator variables may be included in the  unconditional  state- |
      ment.  (See example 3 below.)

 4)   An  IF statement testing ICALL and A_0FLG together is not permit-
      ted.  Instead, two separate nested IF statements must be used: an
      IF statement testing ICALL must occur as the outermost statement,
      and an IF statement testing A_0FLG must occur  as  the  innermost
      statement.   The latter may be supplied by NM-TRAN as a result of
      using an implicit initialization block.  (See example 4 below.)

 5)   Within an explicit compartment initialization block,  A_0(n)  may
      be  assigned conditionally.  The usual rules apply if A_0(n) is a
      random variable. E.g., A_0(n) cannot be assigned within a  nested
      IF,  and  it  defaults  to  0 if it is assigned conditionally but
      incompletely.  However, in checking for a  nested  IF,  tests  of
      A_0FLG and of ICALL are ignored.  (See example 3 below.)

 6)   User-defined  variables may be defined in compartment initializa-
      tion blocks, but not reserved variables such as  basic  or  addi-
      tional PK parameters.

 EXAMPLES OF USAGE:

 The following two fragments of code yield identical results:

 (1) Explicit compartment initialization block:

 IF (A_0FLG.EQ.1) THEN
  A_0(1)=THETA(1)*(1+ETA(1))
 ENDIF

 (2) Implicit compartment initialization block:

 A_0(1)=THETA(1)*(1+ETA(1))

 (3)  This  is  an  example of conditional assignment of A_0(n) (X is a
 data item or user-defined variable):

 IF (A_0FLG.EQ.1) THEN
     IF (X.EQ.1) THEN
     A_0(1)=THETA(1)*(1+ETA(1))
     ELSE
     A_0(1)=THETA(2)*(1+ETA(2))
     ENDIF
 ENDIF

 Note that this can be expressed unconditionally with  the  use  of  an |
 indicator  variable.  E.g.,  if X is a 0/1 variable, then the above is |
 equivalent to                                                          |

     A_0(1)=X*THETA(1)*(1+ETA(1))+(1-X)*A_0(1)=THETA(2)*(1+ETA(2))      |

 (4) Suppose compartment initialization should occur  only  during  the
 Simulation step.  The following is not permitted:

 IF (ICALL.EQ.4.AND.A_0FLG.EQ.1) A_0(1)=THETA(1)*(1+ETA(1))

 Instead, use:

 IF (ICALL.EQ.4) THEN
   IF (A_0FLG.EQ.1) A_0(1)=THETA(1)*(1+ETA(1))
 ENDIF

 or simply

 IF (ICALL.EQ.4) THEN
    A_0(1)=THETA(1)*(1+ETA(1))
 ENDIF

 or even more simply

 IF (ICALL.EQ.4) A_0(1)=THETA(1)*(1+ETA(1))

 (See Compartment Initialization: A_0)
 (See Compartment Initialization: A_0FLG)
 (See pk).

 REFERENCES: None.
