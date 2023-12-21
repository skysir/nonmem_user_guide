


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     YLO YUP PROBABILITY: PR_Y                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_REAL, ONLY: PR_Y

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: PR_Y

 DISCUSSION:

 When with a given data record, either of the limits YLO or YUP are set
 so that during the analysis an interval is defined in which  (or  out-
 side  of  which) an observation is conditioned to exist, then during a
 copying pass (and during ICALL=5 and 6), PR_Y is the estimated  proba-
 bility that an observation will fall within (or outside) the interval.
 (See YLO YUP)
 (See copying_block, expectation)
 (See data average)

 For  the  purpose of computing this probability, an actual observation
 need not exist in the record.  If the mean and variance of the  intra-
 individual model for a potential observation are specified, and if the
 limits YLO or YUP are set, a value of PR_Y will be  computed,  whether
 the  record  has  MDV=0 or MDV=1.  If neither YLO nor YUP are set, the
 value PR_Y will be 1.

 PR_Y may be used as a right-hand  quantity  in  abbreviated  code  for
 $PRED, $PK, $ERROR, and $INFN blocks.

 Location prior to NONMEM 7: rocm44

 REFERENCES: None.
