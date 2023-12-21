


 +--------------------------------------------------------------------+
 |                                                                    |
 |                             MIX: DATA                              |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: MIX routine

 USAGE:
      USE ROCM_REAL, ONLY: DATA=>RDATA

 GLOBAL DECLARATION:
      USE SIZES, ONLY: NO,DPSIZE
      REAL(KIND=DPSIZE) :: RDATA(NO,3)

 DISCUSSION:

  DATA(k,i)
      The  value  of  the  ith type of data item specified in NM-TRAN's
      $CONTR record, appearing on the kth  observation  record  of  the
      individual  record.   It  changes  values  with  each  individual
      record.  DATA is a reserved variable in $MIX abbreviated code.

 DATA is also used in CONTR
 (See CCONTR:_Y,DATA,N1,N2,DIM)

 Location prior to NONMEM 7: rocm1

 REFERENCES: Guide IV, section III.B.4 
