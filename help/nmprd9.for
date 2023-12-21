


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            PASS: PASSRC                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_REAL, ONLY: PASSRC

 GLOBAL DECLARATION:
      USE SIZES, ONLY: VSIZE,DPSIZE
      REAL(KIND=DPSIZE) :: PASSRC(VSIZE)

 DISCUSSION:

  PASSRC
      PASSRC  contains  the  data  record  at  ICALL values 0, 1 and 3.
      PASSRC changes in conjunction with calls to PASS, at  which  time
      transgeneration is permitted.

 Location prior to NONMEM 7: nmprd9

 REFERENCES: None.
