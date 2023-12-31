


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               F_FLAG                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPR_INT, ONLY: F_FLAG=>IPRDFLG1

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IPRDFLG1

 DISCUSSION:

 When  either  of the options LIKELIHOOD or -2LOGLIKIHOOD appear in the
 $ESTIMATION record, this indicates  that  for  all  observations,  the
 variable  Y (with NM-TRAN abbreviated code) or F (with a user-supplied
 PRED or ERROR routine) is being set to a likelihood (or -2 log likeli-
 hood) value for the observation.  Suppose, though, that neither option
 appears.  Then for a given observation, the default indication is that
 Y  or  F  is being set to a "prediction" of the observation.   This is
 equivalent to setting the variable F_FLAG to its default value  of  0.
 When, however, the variable F_FLAG is set to the value 1 (or 2) in the
 PRED or ERROR routine, this signals that Y or F  is  being  set  to  a
 likelihood  (or  -2 log likelihood) value for this particular observa-
 tion.

 F_FLAG may be set in NM-TRAN abbreviated code.

 IF (TYPE.EQ.1) THEN
    Y=THETA(1)+ETA(1)+ERR(1) ; a prediction
 ELSE
    F_FLAG=1
    A=EXP(THETA(2)+ETA(2))
    B=1+A
    Y=DV*A/B+(1-DV)*1/B      ; a likelihood
 ENDIF

 A nonzero value of F_FLAG has no effect during a Simulation Step.  But
 note:  When single-subject data are simulated for which, when the data
 are later analyzed, Y (F) would always be set  to  a  likelihood,  and
 when the ONLYSIMULATION option is used in the $SIMULATION record, then
 usually the NOPREDICTION option can also be used.  When the  NOPREDIC-
 TION  option is used, if any eta's are used in the model, the data are
 regarded as population data, but the model for the  data  in  question
 usually does not involve eta variables.  When, however, single-subject
 data are simulated for which, when the data are analyzed, Y (F)  would
 sometimes  be set to a likelihood and sometimes to a "prediction", and
 when the ONLYSIMULATION option is used in the $SIMULATION record, then
 the  PREDICTION  option  will  need  to  be  used (this is the default
 option!).  The reason is that the expression of  residual  variability
 associated with the predictions will involve eta variables.

 RESTRICTIONS:

 When  F_FLAG  is set to a non-zero-value for some observation, neither
 the LIKELIHOOD or -2LOGLIKIHOOD option may be used.

 When F_FLAG is set to a nonzero-value for an observation,  the  Lapla-
 cian estimation method must be used.

 Unless with all observations within an individual record, F_FLAG is 0,
 the RES and WRES items will be 0 for all data records in the  individ-
 ual record.

 When  F_FLAG  is  set  to a nonzero-value for an observation within an
 individual record, inter-L2 correlated epsilons are not  allowed  with
 the observations within the individual record.
 (See Correlation Across L2 Records)

 If  the  data are population data, a nonzero value of F_FLAG cannot be
 used within a multiple-observation L2 record.  If the data are single-
 subject data, a nonzero value of F_FLAG cannot be used within a multi-
 ple-observation L1 record.

 Location prior to NONMEM 7: nmpr17

 REFERENCES: None.
