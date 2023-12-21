


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        DOSE RECORD: DOSREC                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied PK, DES, AES routines

 USAGE:
      USE PROCM_REAL, ONLY: DOSREC

 GLOBAL DECLARATION:
      USE SIZES, ONLY: VSIZE,DPSIZE
      REAL(KIND=DPSIZE) :: DOSREC(VSIZE)

 DISCUSSION:

 DOSREC
      At  a  call to PK which occurs at a non-event dose time (i.e., at
      the time of an additional or lagged dose),  DOSREC  contains  the
      (last data record of the) event record describing the dose.

      When PK is called at an event time, DOSREC contains 0's.

 (See Dose Time Non-Event: DOSTIM).

 Location prior to NONMEM 7: procm3

 REFERENCES: Guide VI, section III.I 
