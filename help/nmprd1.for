


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           PRED EXIT CODE                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 With versions of NONMEM before 7.4.2:

 USAGE:
      USE NMPRD_INT, ONLY: IERPRD,NETEXT

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IERPRD,NETEXT

 With versions of NONMEM starting with 7.4.2:

 USAGE:
      USE NMPRD_INT, ONLY: IERPRD,IERPRDU,NETEXT

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IERPRD,IERPRDU,NETEXT

 DISCUSSION:

 Values are stored by PRED for use with NONMEM.

 IERPRD
      The  PRED  error return code; also called the PRED exit code.  Is
      set to 0 before PRED is called.
      IERPRD=0: Normal return
      IERPRD=1: PRED is unable to compute.  If possible, NONMEM  should
      attempt recovery.
      IERPRD=2:  PRED  is  unable  to compute.  NONMEM should abort the
      run.

 NETEXT
      NETEXT=0 to 3: the number of lines of text of  an  error  message
      stored by PRED in ETEXT
      (See nmprd2).
      (See PRED Error Message).

 Values may be stored in IERPRD two ways.
 PREDPP may store error messages such as
  PK  PARAMETER  FOR  KA  IS NON-POSITIVE in ETEXT and set IERPRD=1 and
 K=1.  (In effect, EXIT 1 1).

 In abbreviated code, whether or not  PREDPP  is  used,  the  statement
 "EXIT n k" causes n to be stored in IERPRD.  The value of k is part of
 the error message in ETEXT, which is reported  in  the  NONMEM  output
 report and in file PRDERR.

 With NONMEM 7.4.2, there is a new variable, IERPRDU.  Both n and k are
 stored in IERPRDU.  k must be between 0 and 999. The value  stored  is
 n*10000+k.  E.g., "EXIT 1 500" is stored in IERPRDU as 10500.

 With  NONMEM 7.5 and later, values of k between 1000 and 9999 are per-
 mitted.  E.g., "EXIT 1 2000" is stored in IERPRDU as 12000.

 PREDPP will not set these values; only the EXIT statement in  abbrevi-
 ated  code can do this.  IERPRDU is ignored if it is not the Simuation
 Step.  For NONMEM 7.5 use of the error code during simulation,
 (See Simulation Block) and look for "Simulation Error Forgiveness".

 Location prior to NONMEM 7: nmprd1

 REFERENCES: Guide V, section 12.4.15 
 REFERENCES: Guide VI, section III.K.2 , IV.F 
 REFERENCES: Guide IV, section IV.G 
 REFERENCES: Guide Introduction_7
