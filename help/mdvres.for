


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               MDVRES                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY: MDVRES

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: MDVRES

 DISCUSSION:

 MDVRES stands for missing dependent variable (MDV) for residual (RES).
 Setting MDVRES to 1 is equivalent to temporarily declaring an observa-
 tion  as  missing  during  the  computation  of residuals and weighted
 residuals.

 MDVRES=0 (default)

 Set MDVRES to 1 in the $ERROR or $PRED routine if you do not  want  to
 include  a  particular  observation in the computation of residual and
 weighted residuals.

 One situation in which this may be  useful  is  when  F_FLAG  is  used
 because  some  observations  are predictions (F_FLAG=0) and others are
 likelihood or -2 log likelihood values  (F_FLAG=1  or  F_FLAG=2).   By
 default,  if  F_FLAG  is  set  to 1 or 2 for any observation within an
 individual record, the RES and WRES items  will  be  0  for  all  data
 records  in  the   individual record.  When MDVRES is set to 1 for any
 observation, this overrides the default.  NONMEM will compute RES  and
 WRES for all observations for which MDVRES is set to 0.  MDVRES should
 be set to 1 with all observations having F_FLAG=1 or F_FLAG=2, but may
 be set to 1 for other observations as well.

 EXAMPLE:

 Supose  that  some observations are assessed by a non-normal distribu-
 tion likelihood such as the PHI() function for below  detection  limit
 values, in which F_FLAG is set.  By setting MDVRES=1 to these particu-
 lar below detection values, the weighted residual algorithm can assess
 the remaining normally distributed values for that subject.

 $ERROR
 SD = THETA(5)
 IPRED = LOG(F)
 DUM = (LOQ - IPRED) / SD
 CUMD = PHI(DUM)
 IF (TYPE .EQ. 1) THEN
       F_FLAG = 0
       Y = IPRED + SD * ERR(1)
 ENDIF
 IF (TYPE .EQ. 2) THEN
       F_FLAG = 1
       Y = CUMD
       MDVRES=1
 ENDIF

 REFERENCES: Guide Introduction_7
