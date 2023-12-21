


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           WRITE FORMATS                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables

 CONTEXT: User-supplied routines

 USAGE:
      USE PROCM_CHAR, ONLY: FMT

 GLOBAL DECLARATION:
      CHARACTER*80 FMT(3)

 DISCUSSION:

 FMT(1)
      A  suitable  format  statement for WRITE statements when the full
      omega array is written.

 FMT(2)
      A suitable format statement for WRITE statements  when  the  full
      sigma array is written.

 FMT(3)
      A  suitable  format  statement for WRITE statements when the full
      theta array is written.

 Location prior to NONMEM 7: procm8

 REFERENCES: None.
