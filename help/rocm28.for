


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    SUPER PROBLEM PRINT CONTROL                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD_INT, ONLY IPRNV

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IPRNV(2)

 DISCUSSION:

  IPRNV(i)
      0  indicates  no printing of NONMEM input information with itera-
      tions 2, 3, etc. of active ith level superproblem.
      1 indicates printing of NONMEM input information with  iterations
      2, 3, etc. of active ith level superproblem.

 Status with 2nd level superproblem takes precedence.

 (See $super).

 Location prior to NONMEM 7: rocm28

 REFERENCES: None.
