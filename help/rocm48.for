


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     SIZE OF INDIVIDUAL RECORD                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only globoal variable
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_INT, ONLY: LIREC=>NDATPASS

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NDATPASS

 DISCUSSION:
 LIREC is The number of data records contained in the individual record
 at the current call to PRED.

 LIREC also changes value appropriately in conjunction  with  calls  to
 PASS.

 LIREC  may  be used as a right-hand quantity in $PRED, $PK, and $ERROR
 blocks, and in a $INFN block in conjunction with PASS.

 Location prior to NONMEM 7: rocm48

 REFERENCES: None.
