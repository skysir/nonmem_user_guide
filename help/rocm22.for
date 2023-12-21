


 +--------------------------------------------------------------------+
 |                                                                    |
 |                  PARAMETERS OMEGA SIGMA: CURRENT                   |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_REAL, ONLY: OMEGA=>VARNF

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR,DPSIZE
      REAL(KIND=DPSIZE) :: VARNF(LVR,LVR)

 DISCUSSION:

  OMEGA
      The current value of OMEGA being used.

 The current value of SIGMA is also located in this array as follows:
      SIGMA(I,J)=OMEGA(NETAS_+I,NETAS_+J)
 NETAS_ is described in Parameter dimensions.

 (See Parameter dimensions).

 At  run and problem initialization and at problem finalization use the
 variables OMEGAF and SIGMAF.
 (See Parameter values: Initial and Final).

 Should not be used if an initial estimate is being computed.

 When an initial estimate is being computed, the  value  is  just  that
 found
  at problem initialization (when ICALL=1).

 EXAMPLE:

 Compute individual weighted residuals using a slope-intercept residual
 error model:

 $ERROR
 Y=F+EPS(1)+F*EPS(2)
 IF (COMACT.GE.1) THEN
    STD=SQRT(SIGMA(1)+F**2*SIGMA(2))
    IWRES=(DV-F)/STD
 ENDIF

 Location prior to NONMEM 7: rocm22

 REFERENCES: None.
