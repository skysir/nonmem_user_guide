


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                CTLO                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPR_REAL, ONLY:CTLO=>CTLW,DCTLO=>DCTLW,DDCTLO=>DDCTLW

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR,DPSIZE
      REAL(KIND=DPSIZE) :: CTLW,DCTLW(LVR),DDCTLW(LVR,LVR)

 DISCUSSION:
 An  observation may be the event that the value of a normally distrib-
 uted variable falls in a given  interval.   The  likelihood  for  this
 event  may be automatically computed.  With the data record containing
 the observation, one specifies the mean and variance of  the  variable
 for  NONMEM,  as usual, and one sets CTLO to the lower endpoint of the
 interval.  If with population data, this endpoint depends  on  an  eta
 variable,  then the first- and second-derivatives of the endpoint with
 respect to etas should also be set.  Derivatives equal to 0  need  not
 be  explicitly  set,  and the only elements of DDCTLO that need be set
 are those in the lower triangle.  CTLO may be set in a $PRED or $ERROR
 abbreviated code, and then NM-TRAN automatically sets the derivatives.
 NONMEM can detect when CTLO is not set.  When CTLO is set, the  Lapla-
 cian estimation method must be used.

 CTUP may be used in conjunction with CTLO.
 (See ctup).

 PR_CT  is the estimated probability that an observation will be of the
 category in question.
 (See pr_ct, CTLO_CTUP_Probability:_PR_CT).

 Limitation:
 May not be used with the LIKELIHOOD or -2LOGLIKELIHOOD options.

 Location prior to NONMEM 7: nmpr13

 REFERENCES: None.
