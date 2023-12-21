


 +--------------------------------------------------------------------+
 |                                                                    |
 |                 NONPARAMETRIC DENSITY: DEN_,CDEN_                  |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_REAL, ONLY: DEN_=>DEN_NP(1),CDEN_=>DEN_NP(2:LVR+1)

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR,DPSIZE
      REAL(KIND=DPSIZE) :: DEN_NP(LVR+1)

 DISCUSSION:

  DEN_
      The nonparametric density.

  CDEN_(1)
      The marginal cumulative value for the first eta.

  CDEN_(2)
      The marginal cumulative value for the second eta.
      etc.

 Values are computed by NONMEM when the Nonparametric step is performed
 and marginal cumulatives are requested (with  NM-TRAN:  $NONPARAMETRIC
 MARGINALS).

 Values  are  available  during  the pass with COMACT=2 and (if PASS is
 called) at ICALL=3.

 These variables may be used as right-hand  quantities  in  abbreviated
 code blocks $PK, $ERROR, $INFN and $PRED.

 Location prior to NONMEM 7: rocm18

 REFERENCES: None.
