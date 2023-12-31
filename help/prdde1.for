


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     DES AES: ICALL,IDEFD,IDEFA                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: DES and AES routines

 USAGE:
      USE PRMOD_INT, ONLY ICALL=>ICALLD,IDEFD,IDEFA

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: ICALLD,IDEFD(2),IDEFA(2)

 DISCUSSION:

 ICALL
      Identical to the argument ICALL passed by NONMEM to PREDPP.

 IDEFD
      DES may set IDEFD when ICALL is 1, as follows.

      IDEFD(1) may optionally be set by DES to indicate how many thetas
      it uses.  Set to 0 if none.  Otherwise, set to the index  of  the
      highest   numbered  theta  used.   IDEFD(1)=-9  means  "unknown".
      PREDPP determines from the value stored in IDEFD(1) how many ele-
      ments of theta to copy from its input argument THETA to THETAS.
      (See DES_AES:_THETA).

      When  IDEFD(1)  is  set by DES to -9, PREDPP copies all thetas in
      the  problem  to  THETAS  (See DES AES: THETA).   With   NM-TRAN,
      IDEFD(1) is set to -9 when the $DES block contains verbatim code.
      If a user-written DES leaves IDEFD(1) unchanged, it  defaults  to
      -9,  which  does no harm, but may cause the run to be slower than
      necessary.

      IDEFD(2) The full/compact flag for DES.

      DES sets this as follows:
      IDEFD(2)=0 DES will return compact arrays.
      IDEFD(2)=1 DES will return full arrays (the default).

 IDEFA
      AES may set IDEFA when ICALL is 1, as follows.

      IDEFA(1) is set by AES to indicate how many thetas it uses.   See
      above remarks for IDEFD and DES.

      IDEFA(2)  Calling protocol for ADVAN9 and ADVAN15 and ADVAN17

      AES sets this as follows:

      IDEFA(2)=-1 "call with every event record" (the default)
      IDEFA(2)=1 "call once per individual record"

      IDEFA(2) applies only when no TIME data item
      is defined. It is ignored when the TIME data item is defined.

 Location prior to NONMEM 7: prdde1

 REFERENCES: None.
