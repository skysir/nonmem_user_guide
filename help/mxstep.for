


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               MXSTEP                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE PRDATA, ONLY: MXSTEP=>MXSTP0  (ADVAN9)
      USE PRDATA, ONLY: MXSTEP=>MXSTP01 (ADVAN13,ADVAN14,ADVAN15,ADVAN16,ADVAN17,ADVAN18)

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) MXSTP0
      INTEGER(KIND=ISIZE) MXSTP01

 DISCUSSION:

  MXSTEP
      The  maximum  number  of integration steps for ADVAN9 and ADVAN13
      and ADVAN14 and ADVAN15 and ADVAN16 and ADVAN17 and ADVAN18.

      For ADVAN13, ADVAN14, ADVAN15, ADVAN16, ADVAN17, and ADVAN18, the
      default value of MXSTP01 is set to 10000 in resource\PRDATA.f90

      For  ADVAN9,  the  default  value of MXSTP0 is set to the largest
      possible integer  value  2147483647  in  resource\PRDATA.f90,  so
      MXSTEP can only be used to set a smaller value.

      This  variable is only a reserved variable for ADVAN9 and ADVAN13
      and ADVAN14 and ADVAN15 and ADVAN16 and ADVAN17 and ADVAN18.

 REFERENCES: Guide Introduction_7
