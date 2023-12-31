


 +--------------------------------------------------------------------+
 |                                                                    |
 |                PARAMETER VALUES: INITIAL AND FINAL                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_REAL, ONLY: THETAF,OMEGAF,SIGMA,THETAFR

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LTH,LVR ,DPSIZE
      REAL(KIND=DPSIZE) :: THETAF(LTH),OMEGAF(LVR,LVR), &
                           SIGMAF(LVR,LVR),THETAFR(LTH)

 DISCUSSION:

 THETAF
      THETAF(i)  = zero, initial, or final value of theta(i), according
      to the current value of ICALL.

 THETAFR
      THETAFR(i) = zero, initial, or final value of reported thetar(i), |
      according  to  the  current value of ICALL.  If record $THETAR is |
      not used, THETAF and THETAFR are equal.   If  record  $THETAR  is |
      used, then THETAF is the internal theta as used in $PK/$PRED, and |
      THETAFR is the theta reported in the report file.

 OMEGAF
      OMEGAF(i,j) =  zero,  initial,  or  final  value  of  omega(i,j),
      according to the current value of ICALL.

 SIGMAF
      SIGMAF(i,j)  =  zero,  initial,  or  final  value of sigmaf(i,j),
      according to the current value of ICALL.

 At ICALL = 0, these are zero values.

 At ICALL = 1, these are initial values.

 If NONMEM will be computing the initial value of theta(i), at  ICALL=1
 the initial value THETAF(i) is set to some point between the lower and
 upper bounds.
 If NONMEM will be computing the initial value of a diagonal  block  of
 OMEGA (SIGMA), at ICALL=1 the initial value of the block is set to the
 identity matrix.

 At ICALL = 3, these are final values.

 At ICALL = 2, the values of these variables  should not be used unless
 the value of COMACT is 1 or 2, in which case these are final values.
 (See COMACT,COMSAV)

 Location prior to NONMEM 7: rocm6

 REFERENCES: None.
