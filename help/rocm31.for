


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         MIX CONTR: TEMPLT                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied PRED routine

 USAGE:
      USE ROCM_REAL, ONLY: TEMPLT=>VRAW

 GLOBAL DECLARATION:
      USE SIZES, ONLY: VSIZE,DPSIZE
      REAL(KIND=DPSIZE) :: VRAW(VSIZE)

 DISCUSSION:

  TEMPLT
      This variable serves two different purposes.
      At ICALL=6, the template data record is stored in TEMPLT.
      When  MIX  is  called,  the  first  data record of the individual
      record is stored in TEMPLT.

 TEMPLT(n) may be used as a right-hand quantity in blocks  of  abbrevi-
 ated  code  that  test for ICALL=6 and in $MIX blocks.  The "n" may be
 either a position in the data record or the label of  the  data  item,
 e.g. TEMPLT(1) or TEMPLT(ID).

 Location prior to NONMEM 7: rocm31

 REFERENCES: None.
