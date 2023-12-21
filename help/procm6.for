


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           DES AES: THETA                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied DES, AES routines

 USAGE:
      USE PROCM_REAL, ONLY: THETA=>THETAS

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LTH,DPSIZE
      REAL(KIND=DPSIZE) :: THETAS(LTH)

 DISCUSSION:

 THETA
      The  THETA  vector  passed  as an argument by NONMEM to PRED, PK,
      ERROR.
      IDEFD(1) and IDEFA(1) affect how many elements are copied.
      (See DES_AES:_ICALL,IDEFD,IDEFA)

 Location prior to NONMEM 7: procm6

 REFERENCES:  None.
