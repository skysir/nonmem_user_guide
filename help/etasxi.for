


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               ETASXI                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY: ETASXI

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE), ALLOCATABLE, DIMENSION(:)   :: ETASXI

 DISCUSSION:

 With NONMEM 7, an alternative eta shrinkage evaluation using empirical
 Bayes variances (EBVs, or conditional mean  variances)  are  now  also
 reported.   The $ESTIMATION record option ETASTYPE=0 requests that eta
 shrinkage be averaged for all subjects.  This is  the  default.   With
 ETASTYPE=1,  eta  shrinkage  is averaged only among subjects that pro-
 vided a non-zero derivative of their data likelihood with  respect  to
 that eta.

 Reserved  variable  ETASXI(i)  may  be used to specify certain etas of
 particular subjects to be included, or to specify certain etas of cer-
 tain  subjects  to be excluded, from the average eta shrinkage assess-
 ment.  ETASXI stands for eta shrinkage exclude/include.  The subscript
 i refers to ETA(i).

 If ETASXI(i) is set to 2, ETA(i) is included.
 If ETASXI(i) is set to 1, ETA(i) is excluded.

 This overrides whatever would happen based on ETASTYPE.

 EXAMPLE:

 ETASXI(3)=1
 ETASXI(4)=2

 For all subjects, ETA(3) is excluded and ETA(4) is included.

 ETASXI can be set differently for certain subjects.

 EXAMPLE:

 IF(ID==3)  ETASXI(1)=1
 IF(ID==23) ETASXI(3)=2

 This  excludes  ETA(1)  for subject 1, and includes ETA(3) for subject
 23.

 ETASXI may be set only in $PK and $PRED blocks.

 The additional output file root.shm (which stands for  shrinkage  map)
 will contain information which etas were excluded in the eta shrinkage
 assessment.

 REFERENCES: Guide Introduction_7
