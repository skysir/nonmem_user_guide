


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     ADDITIONAL_RECORD_COUNTERS                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      include nonmem_reserved_general

 DISCUSSION:

 In  addition  to  the record counters NDREC and NIREC, NONMEM 7.4 pro-
 vides additional record counters.  All except  IRECIDX   refers  to  a
 record number within the current individual record.

  FIRSTREC
      First record of the individual.

  LASTREC
      Last record of the individual.

  FIRSTOBS
      First  observation  record  of the individual (for which MDV=0 or
      100)

  LASTOBS
      Last observation record of the individual  (for  which  MDV=0  or
      100)

  FIRSTDOS
      With  PREDPP,  first  dose  record  of  the individual (for which
      EVID=1 or 4, and MDV=1 or MDV=101).
      Without PREDPP, -1.

  LASTDOS
      With PREDPP, last dose record of the individual (for which EVID-1
      or 4, and MDV=1 or MDV=101).
      Without PREDPP, -1.

  EFIRSTREC
      First  record  of  the individual during Estimation or Covariance
      Steps (for which MDV=0 or MDV=1)

  ELASTREC
      Last record of the individual  during  Estimation  or  Covariance
      Steps (for which MDV=0 or MDV=1)

  EFIRSTOBS
      First  observation  record of the individual during Estimation or
      Covariance Steps (for which MDV=0)

  ELASTOBS
      Last observation record of the individual  during  Estimation  or
      Covariance Steps (for which MDV=0)

  EFIRSTDOS
      With  PREDPP,  first dose record of the individual during Estima-
      tion or Covariance Steps (for which EVID=1 or 4, and MDV=1).
      Without PREDPP, -1.

  ELASTDOS
      With PREDPP, last dose record of the individual during Estimation
      or Covariance Steps (for which EVID=1 or 4, and MDV=1).
      Without PREDPP, -1.

  IRECIDX
      IRECIDX+1  is  the  starting record number in the NONMEM data set
      for the current individual record.  NDREC starts  at  1  for  the
      first  record of each individual record, so that NDREC+IRECIDX is
      the number of the current data record within the NONMEM data set.

 Like NIREC and NDREC, the additional record counters also change value
 in conjunction with calls to PASS, a NONMEM utility routine.  They may
 all be used as right-hand quantities in $PRED, $PK, and $ERROR blocks,
 and in a $INFN block in conjunction with PASS.

 REFERENCES: Guide Introduction_7
