


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           CONTR: III,DIM                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables

 CONTEXT: User-supplied routines

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 NONMEM  code and previous examples that may be available from advanced
 users.

 USAGE:
      USE ROCM_INT, ONLY:  III=>NL2_RECS, DIM=>L2_DIM

 GLOBAL DECLARATION:
      USE SIZES, ONLY: NO
      INTEGER(KIND=ISIZE) :: NL2_RECS,L2_DIM(NO)

 DISCUSSION:

 These variables change values with each individual record.   They  may
 be used by CONTR.

  III Number of L2 records in the individual record.

  DIM Length of each L2 record.

 When  there  is  no  L2 data item, III = number of observations in the
 individual record, and all values in DIM = 1.

 Location prior to NONMEM 7: rocm2

 REFERENCES: None.
