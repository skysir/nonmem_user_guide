


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         NINDR INDR1 INDR2                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_INT,  ONLY: NINDR=>NINDOBS,INDR1=>IDXOBSF,INDR2=>IDXOBSL

 GLOBAL DECLARATION:
      REAL(KIND=ISIZE):: NINDOBS,IDXOBSF,IDXOBSL

 DISCUSSION:

  NINDR
      The  number  of  individual records in the data set containing an
      observation record.

  INDR1
      The index of the first individual record in the data set contain-
      ing an observation record.

  INDR2
      The  index of the last individual record in the data set contain-
      ing an observation record.

 These variables may be used as right-hand quantities  in  $PRED,  $PK,
 $ERROR, $INFN and $MIX abbreviated code.

 Location prior to NONMEM 7: rocm46

 REFERENCES: None.
