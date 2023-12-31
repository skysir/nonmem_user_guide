


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               ICALL                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Argument of user-supplied subroutines
 CONTEXT: User-supplied routines

 DISCUSSION:
 ICALL  is  an  argument  passed by NONMEM to user-supplied subroutines
 CCONTR, CONTR, CRIT, PRED, PRIOR and MIX.

 PREDPP selectively passes the same argument to PK,  ERROR,  and  INFN.
 ICALL  can  be  used  in $PK, $ERROR, and $PRED abbreviated code.  The
 discussion below describes the values of ICALL as seen by PRED.

 ICALL=-1: the routine has been called for the PRED_IGNORE_DATA feature
 of  NONMEM  7.5.   One  call per data record, at the start of the run.
 These  calls  occur  only   if   abbreviated   code   uses   variables
 PRED_IGNORE_DATA or PRED_IGNORE_DATA _TEST, or if the PRED_IGNORE_DATA
 option of $DATA is used.  Otherwise, there are no calls to  PRED  with
 ICALL=-1.

 ICALL=0:  the routine has been called for initialization at the begin-
 ning of the NONMEM run; one such call per run.

 ICALL=1: the routine has been called for initialization at the  begin-
 ning of a NONMEM problem; one such call per problem.

 ICALL=2:  the routine has been called for a prediction. Multiple calls
 occur.

 ICALL=3: the routine has been called for finalization at the end of  a
 NONMEM problem; one such call per problem.

 ICALL=4:  the routine has been called during the Simulation Step; mul-
 tiple calls occur, as with ICALL=2.

 ICALL=5: the routine has been called when expectations are being  com-
 puted; multiple calls occur.

 ICALL=6:  the routine has been called when raw data averages are being
 computed; multiple calls occur.

 Note: Some subroutines are called with only a subset of  the  possible
 values of ICALL.

 REFERENCES: Guide  I, section C,4.2 
 REFERENCES: Guide IV, section IV.D , IV.E.1 
 REFERENCES: Guide VI, section III.C , IV.B.1 
