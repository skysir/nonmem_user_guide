


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      MODEL EVENT TIME: MTDIFF                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: PK routine

 USAGE:
      USE PRMOD_INT, ONLY: MTDIFF

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: MTDIFF

 DISCUSSION:

 MTDIFF is of interest when model time parameters MTIME are used.

 The  value  of  MTDIFF is 0 when PK is called.  If PK sets MTDIFF to a
 value other than 0, e.g., MTDIFF=1, then PREDPP will  understand  that
 with  that  call to PK, the values of one or more of the MTIME(i) have
 possibly been reset.

 It is not necessary to set MTDIFF at a  call  to  PK  with  the  first
 record  of an individual or with a reset record.  At such calls, it is
 assumed that all MTIME(i) are being set  (if  only  to  their  default
 value of 0).

 Suppose PK has defined a value MTIME(i).  PREDPP calls PK prior to the
 advance to MTIME(i).  PK may change MTIME(i) at this  call.   However,
 the  new  value  of  MTIME(i) is ignored till after the advance.  Only
 then is the new value of MTIME(i) effective.

 MTDIFF=0 (the default) can save considerable run time when  there  are
 many  model  time parameters.  Note that the results are unpredictable
 if the times are in fact changed when MTDIFF=0.

 (See mtime, model time examples).

 Location prior to NONMEM 7: prdpk2

 REFERENCES: None.
