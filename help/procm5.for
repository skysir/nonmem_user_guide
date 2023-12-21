


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     ACTIVE ETA LIST FOR PREDPP                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: PK, ERROR, TRANS routines

 USAGE:
      USE PROCM_INT, ONLY: NACTIV,M=>IDXETA

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR
      INTEGER(KIND=ISIZE) :: NACTIV,IDXETA(LVR)

 DISCUSSION:

 NACTIV
      NACTIV = # of etas in the problem - NRETA.
      (See Non-Active_ETA_list_for_PRED).

      NACTIV tells how many of the etas in the problem have partial de-
      rivatives that NONMEM is not currently ignoring.

 M    M(k) is the index of the  kth  0-valued  element  of  LVOUT  (for
      k=1,...,NACTIV).

 Suppose  a user-written TRANS routine modifies the kth. element of GG.
 Here is code that might be used to modify only those of its first  eta
 derivatives that are of current interest to NONMEM:

         DO 100 I=1,NACTIV
     100 GG(K,M(I)+1,1)=GG(K,M(I)+1,1)  + ...

 Location prior to NONMEM 7: procm5

 REFERENCES:  None.
