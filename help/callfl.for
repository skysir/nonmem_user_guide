


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               CALLFL                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Control calls to PK and ERROR subroutine and ADVAN9/15/17
 CONTEXT: Abbreviated code

 USAGE:
 CALLFL=n

 DISCUSSION:
 This  pseudo-assignment  statement  may  be present in $PK abbreviated
 code, in which case it controls when PREDPP calls the PK routine.   It
 may  be  present in $ERROR abbreviated code, in which case it controls
 when PREDPP calls the ERROR routine.  It may be present in $AESINITIAL
 abbreviated code when the TIME data item is not used, in which case it
 controls when PREDPP calls the ADVAN9, ADVAN15,  or  ADVAN17  routine.
 The use of CALLFL in each of these three abbreviated codes is indepen-
 dent of its use in the others.  The pseudo-assignment statement may be
 enclosed in parentheses.  A calling protocol phrase may be used within
 parentheses instead of a pseudo-assignment statement CALLFL (See call-
 ing protocol).  The pseudo-assignment statement may take these forms:

 CALLFL=-2

  Call  the  PK subroutine with every event record, with additional and
  lagged doses, and at modeled event times.   Does  not  apply  to  the
  ERROR or ADVAN9, ADVAN15, or ADVAN17 routines.

  The following messages will appear in the NONMEM output report:
  PK SUBROUTINE CALLED WITH EVERY EVENT RECORD
  PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE TIMES
  or
  PK SUBROUTINE CALLED AT NONEVENT (ADDITIONAL AND LAGGED) DOSE TIMES
  AND AT MODEL TIMES

 CALLFL=-1

  Call the subroutine with every event record.  This is the default.

  Some of the following messages will appear in the NONMEM output:
  PK SUBROUTINE CALLED WITH EVERY EVENT RECORD
  PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES

  ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD

  ADVAN9 CALLED WITH EVERY EVENT RECORD
  ADVAN15 CALLED WITH EVERY EVENT RECORD
  ADVAN17 CALLED WITH EVERY EVENT RECORD

 CALLFL=0

  For the PK subroutine: If the data are population data, call the sub-
  routine with the first event record of each individual record; if the
  data  are  single-subject  data,  call  the subroutine with the first
  event record of the data set.  In addition, call the subroutine  with
  with  every event record where the event time differs from the previ-
  ous event time.

  The following messages will appear in the NONMEM output report:
  PK SUBROUTINE CALLED ONLY WITH NEW INDIVIDUAL OR NEW TIME
  PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL DOSE OR LAGGED) DOSE TIMES.

  For the ERROR subroutine: When the Simulation Step  is  being  imple-
  mented, call the subroutine with every event record.  Otherwise, call
  the subroutine only with observation event records.

  The following messages will appear in the NONMEM output report:
  DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD
  OTHERWISE, ERROR SUBROUTINE CALLED ONLY WITH OBSERVATION EVENTS

  Does not apply to the ADVAN9, ADVAN15, or ADVAN17 routine.

  With CALLFL=0, The CALL data item may be used to request  calls  with
  additional event records.

 CALLFL=1

  For  the  PK and ADVAN9, ADVAN15, or ADVAN17 subroutines: If the data
  are population data, call the subroutine with the first event  record
  of  each individual record; if the data are single-subject data, call
  the subroutine with the first event record of the data set.

  The following messages appear in the NONMEM output report:
  PK SUBROUTINE CALLED ONCE PER INDIVIDUAL RECORD
  PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL DOSE OR LAGGED) DOSE TIMES

  ADVAN9 CALLED ONCE PER INDIVIDUAL RECORD.
  ADVAN15 CALLED ONCE PER INDIVIDUAL RECORD.
  ADVAN17 CALLED ONCE PER INDIVIDUAL RECORD.

  For the ERROR subroutine: When the Simulation Step  is  being  imple-
  mented,  call  the subroutine with every event record.  Otherwise, if
  the data are population data, call  the  subroutine  with  the  first
  event  record  of each individual record; if the data are single-sub-
  ject data, call the subroutine with the first  event  record  of  the
  data set.

  The following messages appear in the NONMEM output report:
  DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD
  OTHERWISE, ERROR SUBROUTINE CALLED ONCE PER INDIVIDUAL RECORD

 With  CALLFL=1,  The  CALL data item may be used to request calls with
 additional event records.

 In a block of abbreviated code, the CALLFL=n pseudo-assignment  state-
 ment  must precede all the other abbreviated code (except for verbatim
 code or other pseudo-assignment  statements).   The  pseudo-assignment
 statement  may not be used conditionally.  CALLFL may not be used as a
 variable elsewhere in the abbreviated code.

 Note: If the $ERROR record consists  of  exactly  one  of  these  four
 statements:
  Y=F+ERR(1)
  Y=F*(1+ERR(1))
  Y=F+F*ERR(1)
  Y=F*EXP(ERR(1))
 with  no  other lines of code (except for comment lines), NM-TRAN will
 automatically limit calls to ERROR  to  once-per-problem  (unless  the
 Simulation  Step  is being implemented, in which case.  calls are made
 with every event record).  In effect, this amounts to yet another  way
 to  control  when it is that calls may occur to the ERROR routine, but
 one which may not be explicitly specified in $ERROR  via  the  use  of
 CALLFL.

 With  the  last  three  models (proportional and exponential), NM-TRAN
 will also cause PREDPP to output the message:
 ERROR IN LOG Y IS MODELED
 (This does not mean that a model is fit to Log Y data.)

 REFERENCES: Guide IV, section V.C.5 
