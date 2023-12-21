


 +--------------------------------------------------------------------+
 |                                                                    |
 |                   SIMULATION: SIMEPS ERROR CODE                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine
 USAGE:
      USE NMPRD_INT, ONLY: IERSQ

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IERSQ

 DISCUSSION:

 Value  is  stored  by SIMEPS for use with PRED and ERROR.  Is relevant
 when correlations are stored in CORRL2.  Is set to 0 before SIMEPS  is
 called.
 (See Correlation_Across_L2_Records).

 IERSQ
      The SIMEPS error return code.
      IERSQ=0: Normal return
      IERSQ=1:  Correlation  matrix  obtained  using  the  correlations
      stored in CORRL2 is not positive definite.

      PRED may attempt corrective action; if it does  so  successfully,
      it should reset IERSQ to 0.

 (With versions of NONMEM through 7.2, C was used rather than CORRL2.)  |

 With  NONMEM 7.3 (and all earlier versions), the value of IERSQ is not |
 checked in generated code.

 Location prior to NONMEM 7: nmprd6

 REFERENCES: None.
