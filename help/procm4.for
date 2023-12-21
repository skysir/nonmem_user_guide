


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          STATE VECTOR: A                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied PK and ERROR routines

 USAGE:
      USE PROCM_REAL, ONLY A=>AMNT,DAETA,D2AETA

 GLOBAL DECLARATION:
      USE SIZES, ONLY: PC,LVR,DPSIZE
      REAL(KIND=DPSIZE) :: AMNT(PC),DAETA(PC,LVR),D2AETA(PC,LVR,LVR)

 DISCUSSION:

 A    A is the state vector of compartment amounts.
      A(n) = the amount in compartment n.

 DAETA
      DAETA(n,i) = the derivative of A(n) wrt eta(i).

 D2AETA
      D2AETA(n,i,j)  = the second derivative of A(n) wrt eta(i), eta(j)
      (lower-triangular; j=1, ..., i).

 The A(n) can be used as right-hand quantities in $ERROR and $PK abbre-
 viated  code.    These  amounts  are  the  latest ones computed.  With
 $ERROR, the A(n) are computed at the event time on  the  event  record
 passed to ERROR.  With $PK, the A(n) are computed at the event time on
 the previous event record, or possibly at a later time.  This  time  -
 the  latest  time  at which the amounts are computed - is given in the
 variable TSTATE, which may also be used in $PK abbreviated code.
 (See State_Vector_TIME:_TSTATE).

 Location prior to NONMEM 7: procm4

 REFERENCES: Guide VI, section IV.D , Figure 14
