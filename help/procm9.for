


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     STATE VECTOR TIME: TSTATE                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied PK and ERROR routine

 USAGE:
      USE PROCM_REAL, ONLY: TSTATE

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: TSTATE

 DISCUSSION:

 TSTATE
      TSTATE = the time at which the state-vector (See State Vector: A)
      was last computed.  It is the previous event time (i.e. the  time
      on  the  previous  event  record  passed to PK), or if at a later
      time, but before the time for which PK is being called, a  lagged
      or  additional  dose  was given, or a regular infusion was termi-
      nated, or a modeled event occurred, then  TSTATE  is  the  latest
      such time.

 NM-TRAN  includes  this  global  variable in the PK and ERROR routines
 when $PK or $ERROR abbreviated code includes references  to  variables
 A(n) or when verbatim code is present.

 Location prior to NONMEM 7: procm9

 REFERENCES: None.
