


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     PASS NEW L2 RECORD: NEWL2                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD_INT, ONLY: NEWL2

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NEWL2

 DISCUSSION:

 NEWL2
      NEWL2 = 1 if the data record is the first of an L2 record.
      NEWL2 = 2 otherwise.

 If  there  is  no  L2  data  item,  each data record is a different L2
 record.

 NEWL2 also changes value in conjunction with calls to PASS.

 NEWL2 may be used as a right-hand quantity in $PRED and $ERROR blocks,
 and in an $INFN block in conjunction with PASS.

 Location prior to NONMEM 7: rocm17

 REFERENCES: None.
