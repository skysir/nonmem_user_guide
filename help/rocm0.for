


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          MIX CONTR: THETA                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables

 CONTEXT: MIX, CONTR, CCONTR routine

 USAGE:
      USE ROCM_REAL, ONLY: THETA=>THETAC

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LTH,DPSIZE
      REAL(KIND=DPSIZE) :: THETAC(LTH)

 DISCUSSION:

  THETA
      The current value of theta, made available by NONMEM for the MIX,
      CONTR and CCONTR routines.  THETA is a reserved varible  in  $MIX
      abbreviated code.

 Location prior to NONMEM 7: rocm0

 REFERENCES: Guide VI, section III.L.2 , Figure 6
