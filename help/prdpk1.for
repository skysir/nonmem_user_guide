


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      MODEL EVENT TIME: MTIME                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: PK routine

 USAGE:
      USE PKERR_REAL, ONLY :: MTIME

 GLOBAL DECLARATION:
      USE SIZES, ONLY: PCT,DPSIZE
      REAL(KIND=DPSIZE) :: MTIME(PCT)

 DISCUSSION:

 Model time parameters MTIME are stored in MTIME if they are defined in
 $PK abbreviated code or if verbatim code  is  present  in  $PK.   This
 makes  them  also  available to the ERROR routine.  However, their eta
 derivatives (if any) are not included.  Hence the ERROR routine should
 not  use model time parameters in such a way as to influence the value
 of Y if those parameters have eta derivatives.

 A user-written PK routine need not define these variables unless  per-
 haps the model time parameter values are used in the ERROR routine.

 (See mtime, model time examples).

 Location prior to NONMEM 7: prdpk1

 REFERENCES: None.
