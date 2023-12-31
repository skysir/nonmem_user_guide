


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $CONTR                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Defines values for certain user-supplied routines
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $CONTR   DATA=([label1|0] [label2|0] [label3|0])

 The  data  item with the Jth label (J=1,2,3) and from the Ith observa-
 tion record of an individual record is available in DATA(I,J).   If  0
 is used instead of a label, then a zero appears in DATA(I,J).

 SAMPLE:
 $CONTR    DATA=(0,TYPE)

 DISCUSSION:
 Optional.  Used only with user-supplied routines such as MIX and CONTR
 and CCONTR that use data items stored in the DATA array.  This  record
 gives  labels  (or  synonyms)  defined  in the $INPUT record of one to
 three types of data items to be made available to the subroutine(s) in
 the DATA array.  These routines are called with individual records. An
 array DATA is available in NONMEM module ROCM_REAL and  changes  value
 with each individual record.

 With  the  above  sample  $CONTR  record,  the following code might be
 present in a double precision MIX routine.  The code loops through the
 observation records of the NREC'th individual record.  For each of the
 NOBS observation records, the local variable TYPE is given  the  value
 of the TYPE data item from that data record.  The 0 in the sample is a
 place-holder which causes the first column in the  DATA  array  to  be
 skipped.   The  value of TYPE for the Ith observation record is there-
 fore available in DATA(I,2).  The DATA array is  found  in  ROCM_REAL.
 NO  is a constant giving the maximum number of observations per  indi-
 vidual  record.  NOBS is the number of  observations  in  the  current
 individual record.  (See sizes).

      USE SIZES, ONLY: NO,DPSIZE
      USE ROCM_REAL, ONLY: DATA=>RDATA
      USE ROCM_INT, ONLY: NOBS=>NOBSIND2
      ...
      INTEGER I
      REAL(KIND=DPSIZE) :: TYPE
      ...
      DO 100 I=1,NOBS
      TYPE=DATA(I,2)
      ...
  100 CONTINUE
      ...

 REFERENCES: Guide IV, section III.B.4 
