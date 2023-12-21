


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        MIXTURE MODEL: MIXPT                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_REAL, ONLY: MIXPT=>MIXP_RAW

 GLOBAL DECLARATION:
      USE SIZES, ONLY: MMX,DPSIZE
      REAL(KIND=DPSIZE) :: MIXP_RAW(MMX)

 DISCUSSION:

  MIXPT
      At  ICALL=6,  (the  final  estimate of) the mixture probabilities
      associated with the individual  record  containing  the  template
      data record are stored in MIXPT.

 Location prior to NONMEM 7: rocm47

 REFERENCES: None.
