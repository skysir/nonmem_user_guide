


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      CCONTR: Y,DATA,N1,N2,DIM                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 NONMEM  code and previous examples that may be available from advanced
 users.

 USAGE:
      USE ROCM_REAL, ONLY:  Y=>DV_ITM2,DATA=>CDATA
      USE ROCM_INT, ONLY: N1=>NL2_OBS,N2=>L2NO,DIM=>L2CR_DIM

 GLOBAL DECLARATION:
      USE SIZES, ONLY: NO,DPSIZE
      REAL(KIND=DPSIZE) :: DV_ITM2(NO),CDATA(NO,3)
      INTEGER(KIND=ISIZE) :: NL2_OBS,L2NO,L2CR_DIM(NO)

 DISCUSSION:

 These variables change values with each L2 record.  They may  be  used
 by CCONTR.
 When  the  epsilons  are correlated between L2 records these variables
 change values with each L1 record.
 (See Correlation_Across_L2_Records)

  Y(k)
      DV data item on the kth observation record of the L2 (L1) record.

  DATA(k,i)
      The ith type of data item specified in NM-TRAN's  $CONTR  record,
      appearing on the kth observation record of the L2 (L1) record.

  N1  Number  of  observations in the L2 record.  When the epsilons are
      correlated between L2 records, N1 is the  number  of  L2  records
      within the L1 record.

  N2  Number  of the L2 record within the L1 record.  When the epsilons
      are correlated between L2 records, N2 is 0.

  DIM When the epsilons are correlated between L2  records,  DIM  gives
      the lengths of the L2 records.

 Location prior to NONMEM  7: rocm3

 REFERENCES: None.
