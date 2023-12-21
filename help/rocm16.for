


 +--------------------------------------------------------------------+
 |                                                                    |
 |                 SIGNIFICANT DIGITS FROM EST. STEP                  |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_REAL, ONLY NSIG=>SIGD,SIG=>DIFA

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LPAR,DPSIZE
      REAL(KIND=DPSIZE) :: SIGD,DIFA(LPAR)

 DISCUSSION:

  SIG(i)

 When  the minimization in the Estimation Step terminates due to round-
 ing errors, but the number of significant digits that is  achieved  is
 reported,  the ith value in SIG gives the number of significant digits
 for the ith UCP element.

  NSIG

 NSIG gives the minimum of the values found in the vector SIG

 The NSIG and SIG values should only be used with ICALL = 3.

 Location prior to NONMEM 7: rocm16

 REFERENCES: None.
