


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      OBJECTIVE FUNCTION VALUE                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_REAL, ONLY: OBJECT

 GLOBAL DECLARATION
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: OBJECT

 DISCUSSION:

 OBJECT
      The final value of the objective function.

 This value should only be used at ICALL = 3 (finalization block).

 May be used as a right-hand quantity in abbreviated code.

 Location prior to NONMEM 7: rocm8

 REFERENCES: None.
