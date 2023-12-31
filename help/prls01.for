


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      STOP_TIME,ITASK_ (NM75)                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE PRLS01_REAL, ONLY: STOP_TIME
      USE PRLS01_INT, ONLY: ITASK_=>ITASK

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL (KIND=DPSIZE) STOP_TIME
      USE SIZES, ONLY: ISIZE
      INTEGER(KIND=ISIZE) ITASK_

 DISCUSSION:

 These  reserved  variables  are  used  to  avoid  overshoot in ADVAN9,
 ADVAN13, ADVAN14, and ADVAN15.  These  LSODA  based  routines  use  an
 algorithm  for  integration  that  overshoots the integration interval
 during calls to DES, but still accurately evaluates at the end of  the
 integration interval when all calculations are completed. However, you
 may wish to capture a maximal or minimal value during  $DES,  and  the
 overshoot  should be turned off for this purpose. This is readily done
 by setting ITASK_=4 in $PK or $ERROR. E.g.,

 $PK ITASK_=4

 ITASK_ may take values between 1 and 5.

 For other values of ITASK_ and a discussion of these variables,
 See INTRODUCTION TO NONMEM 7, "ITASK_ and  STOP_TIME:  Avoiding  over-
 shoot in ADVAN9, ADVAN13, ADVAN14, and ADVAN15"

 You  may  also  specify  a  STOP_TIME (Tcrit) past which it should not
 integrate, if it is different from the end of the  normal  integration
 interval:

 IF(TIME==4.0) STOP_TIME=5.0

 To set back to default (end of normal integration interval),

 STOP_TIME=-1.0d+300

 REFERENCES: Guide Introduction_7
