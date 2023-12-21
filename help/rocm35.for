


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        PARAMETER DIMENSIONS                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD_INT, ONLY: NTHES_=>NWTHT,NETAS_=>NWETA, NEPSS_=>NWEPS

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NWTHT,NWETA,NWEPS

 DISCUSSION:

  NTHES_
      The dimension of theta. If the dimension is 0, NTHES_=1.

  NETAS_
      The dimension of omega. If the dimension is 0, NETAS_=1.

  NEPSS_
      The dimension of sigma. If the dimension is 0, NEPSS_=1.

 Location prior to NONMEM 7: rocm35

 REFERENCES: None.
