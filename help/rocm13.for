


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            CONTR: KCALL                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 NONMEM  code and previous examples that may be available from advanced
 users.

 USAGE:
      USE ROCM_INT, ONLY: KCALL=>K_CONTR

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: K_CONTR

 DISCUSSION:

  KCALL
      KCALL=1 with a unique pass through data calling  PRED  and  CONTR
      with final parameter estimates.
      KCALL=0 otherwise.

 Location prior to NONMEM 7: rocm13

 REFERENCES: None.
