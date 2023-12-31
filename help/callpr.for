


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      CALLING PROTOCOL PHRASE                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING:  Control  calls  to  PK  and  ERROR  and ADVAN9, ADVAN15, and
 ADVAN17
 CONTEXT: Abbreviated code

 USAGE:
 $PK (ONCE PER INDIVIDUAL RECORD)
 $ERROR (ONCE PER INDIVIDUAL RECORD)
 $AESINITIAL (ONCE PER INDIVIDUAL RECORD)

 DISCUSSION:
 A calling protocol phrase can be used instead of  the  CALLFL  pseudo-
 statement  in  $PK,  $ERROR,  and  $AESINITIAL abbreviated codes.  The
 phrase must be enclosed in parentheses.  Either upper  or  lower  case
 may  be used.  In an abbreviated code, the line of code containing the
 phrase must precede all the other abbreviated code (except for  verba-
 tim  code  or  other  pseudo-assignment statements), and it may be the
 same line that marks the beginning of the code, as in the above  usage
 examples.   Pseudo-statements  defining COMRES may be coded within the
 same parentheses, separated by a semicolon ";".  No  abbreviated  code
 may follow ")" on the same line as the ")".

 Phrases equivalent to CALLFL=-2 for the $PK record:

  (NON-EVENT)
  (ADDITIONAL OR LAGGED)

  The following messages will appear in the NONMEM output report:

  PK SUBROUTINE CALLED WITH EVERY EVENT RECORD
  PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE TIMES
  or
  PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE
  TIMES AND AT MODEL TIMES

 Phrases  equivalent  to  CALLFL=-1 for the $PK, $ERROR and $AESINITIAL
 records (the default):

  (EVERY EVENT)
  (EVERY)

  The following messages will appear in the NONMEM output report (as appropriate):

  PK SUBROUTINE CALLED WITH EVERY EVENT RECORD
  PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES

  ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD

  ADVAN9 CALLED WITH EVERY EVENT RECORD.
  ADVAN15 CALLED WITH EVERY EVENT RECORD.

 Phrases equivalent to CALLFL=0 for the $PK record:

  (NEW TIME)
  (NEW EVENT TIME)

  The following messages appear in the NONMEM output report:

  PK SUBROUTINE CALLED ONLY WITH NEW INDIVIDUAL OR NEW TIME
  PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.

 Phrases equivalent to CALLFL=0 for the $ERROR record:

  (OBSERVATION EVENT)
  (OBS)
  (OBSERVATION ONLY)
  (OBS ONLY)

  One or more of the following messages appear  in  the  NONMEM  output
  report when calls are limited in this manner:

  DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD
  OTHERWISE, ERROR SUBROUTINE CALLED ONLY WITH OBSERVATION EVENTS

 Phrases  equivalent  to  CALLFL=1  for the $PK, $ERROR and $AESINITIAL
 blocks:

  (ONCE PER INDIVIDUAL RECORD)
  (ONCE/IND.REC.)

  The following messages appear in the NONMEM output (as appropriate):

  PK SUBROUTINE CALLED ONCE PER INDIVIDUAL RECORD
  PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES

  DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD
  OTHERWISE, ERROR SUBROUTINE CALLED ONCE PER INDIVIDUAL RECORD

  ADVAN9 CALLED ONCE PER INDIVIDUAL RECORD.
  ADVAN15 CALLED ONCE PER INDIVIDUAL RECORD.

 REFERENCES: Guide IV, section V.C.5 
