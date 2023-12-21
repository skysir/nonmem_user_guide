


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      ESTIM COVAR ERROR CODES                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables

 CONTEXT: PRED routine

 USAGE:
      USE ROCM_INT, ONLY: IERE=>IEST_ERR,IERC=>ICOV_ERR

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IEST_ERR,ICOV_ERR

 DISCUSSION:

 IERE The return code from the Estimation Step.

 IERC The return code from the Covariance Step.

 Values of 0 indicate normal termination.

 These  variables  may  be used as right-hand quantities in abbreviated
 code for initialization/finalization blocks.

 Location prior to NONMEM 7: rocm9

 REFERENCES: None.
