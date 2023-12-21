


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    RECORD COUNTERS: NIREC,NDREC                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_INT, ONLY: NIREC=>NINDREC,NDREC=>NDATINDR

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NINDREC,NDATINDR

 DISCUSSION:

  NIREC
      The number of the individual record at the current call.

  NDREC
      The number of the data record within the individual record at the
      current call.

 NIREC and NDREC also change value in conjunction with calls to PASS, a
 NONMEM utility routine.

 NIREC  and  NDREC  may be used as right-hand quantities in $PRED, $PK,
 and $ERROR blocks, and in a $INFN block in conjunction with PASS.

 NM-TRAN automatically provides the necessary USE statement in the gen-
 erated subroutines.

 See  also  additional_record_counters  for other variables of interest
 such as FIRSTREC.

 Location prior to NONMEM 7: rocm32

 REFERENCES: None.
