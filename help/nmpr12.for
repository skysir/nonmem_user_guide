


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              YLO YUP                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPR_REAL, ONLY: YLO,YUP

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: YLO,YUP

 DISCUSSION:
 If an observation is assumed to be normally distributed, with mean and
 variance specified for NONMEM as usual, then the  likelihood  for  the
 observation is conditioned on the observation being in (or outside) an
 interval defined by the values YLO and YUP.  If YLO is less than  YUP,
 then  the likelihood of the observation is conditioned on the observa-
 tion being in the interval with lower bound YLO and upper  bound  YUP.
 If  YLO is greater than YUP, then the likelihood of the observation is
 conditioned on the observation being outside the interval  with  lower
 bound  YUP and upper bound YLO.  These values may be set with the data
 record containing the observation.  They  may  be  set  in  $PRED  and
 $ERROR  abbreviated codes.  If YLO or YUP is not set, it is assumed to
 be minus infinity or plus infinity, respectively.  When YLO or YUP  is
 set, the Laplacian estimation method must be used.

 PR_Y is the estimated probability that an observation will fall within
 (or  outside) the  interval.

 (See pr_y, YLO_YUP_Probability:_PR_Y).

 Location prior to NONMEM 7: nmpr12

 REFERENCES: None.
