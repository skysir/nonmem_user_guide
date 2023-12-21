


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          DES AES: ISFINL                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied DES and AES routines

 USAGE:
      USE PROCM_INT, ONLY: ISFINL

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: ISFINL

 DISCUSSION:

 ISFINL
      When   NONMEM   is   performing  simulation  or  a  copying  pass
      (COMACT>0), DES and AES are called immediately after the  advance
      to  an  event  time or non-event time, with a value of T equal to
      this time.

      ISFINL = 1 at such a call; otherwise =0.

 NM-TRAN includes this global variable in the DES or AES  routine  when
 the $DES or $AES block includes references to variable ISFINL, or when
 verbatim code is present.

 Location prior to NONMEM 7: procmb

 REFERENCES: None.
