


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         PASS NEWIND: NWIND                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD_INT, ONLY: NWIND

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NWIND

 DISCUSSION:

  NWIND
      The  NEWIND  value  at  ICALL 0, 1 and 3.  NWIND changes value in
      conjunction with calls to PASS.  (See pred, newind).

 Location prior to NONMEM 7: rocm34

 REFERENCES: Guide I, section C.3.5.2 
