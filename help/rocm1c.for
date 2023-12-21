


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         CONTR: Y DATA NOBS                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: CONTR routine

 USAGE:
      USE ROCM_REAL, ONLY: Y=>DV_ITM,DATA=>RDATA
      USE ROCM_INT,  ONLY: NOBS=>NOBSIND2

 GLOBAL DECLARATION:
      USE SIZES, ONLY: NO,DPSIZE
      REAL(KIND=DPSIZE) :: DV_ITM(NO),RDATA(NO,3)
      INTEGER(KIND=ISIZE) :: NOBSIND2

 DISCUSSION:

 These  variables  change values with each individual record.  They may
 be used by a user-written CONTR subroutine.

  Y(k)
      DV data item on the kth data record  of  the  individual  record,
      ignoring data records with MDV=1.

  DATA(k,i)
      The  value  of  the  ith type of data item specified in NM-TRAN's
      $CONTR record, appearing on the kth  observation  record  of  the
      individual record.

  NOBS
      Number of observations in the individual record.

 DATA is used also in MIX (See data_for_mix).

 Location prior to NONMEM 7: rocm1

 REFERENCES: Guide IV, section III.B.4 
