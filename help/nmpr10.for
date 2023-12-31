


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        REPETITION VARIABLES                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPR_INT, ONLY: RPTI=>NRPT_IN,RPTO=>NRPT_OUT, &
                          RPTON=>NRPT_ON,PRDFL=>IUSEPRD

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NRPT_IN,NRPT_OUT,NRPT_ON,IUSEPRD

 DISCUSSION:

 These  variables  provide  the  information controlling the Repetition
 feature of NONMEM.

  RPTO
      RPTO is the "repetition output value".  With each data record, it
      may  be  set by PRED, and it conveys information to NONMEM.  RPTO
      may be used as a left-hand quantity  in  $PRED,  $PK  and  $ERROR
      blocks.

       RPTO=0
           If RPTO is set to 0, this produces no effect.

       RPTO=n (where n is between 1 and 5)
           The current record is thus marked as a "repetition base with
           value n", i.e.  as the  first  of  a  series  of  contiguous
           records  of  the current individual record (with single-sub-
           ject data, contiguous records of the data set) which will be
           repeated.  See the repeat data item (See RPT_).

       RPTO=-n (where n is between 1 and 5)
           The current record thus becomes a "repetition initiator with
           value n", i.e.  before the next record is passed to PRED all
           preceding records of the individual record (with single-sub-
           ject data, all records of the data set), starting  with  the
           last  such  record marked as a repetition base with value n,
           and up to and including the current  record,  will  be  once
           again  passed  to PRED (will be "repeated").  The value n is
           put onto the top of a stack - the "repetition stack".

           The series of records from base to initiator is  called  the
           "repetition  series  for n".  The series may be repeated any
           number of times (see RPTON).  After all repetitions are com-
           plete, the value n is removed from the repetition stack.

           While  the  series  is being repeated, RPTO may be set again
           with any record R of the series.  If the (unsigned)  repeti-
           tion  output  value  is  a  value  already on the repetition
           stack, this value is ignored.  Otherwise, the  output  value
           is not ignored, and it may be set so as to mark R as a repe-
           tition base or as a repetition  initiator.   In  the  latter
           case the current value on the top of the stack (m) is pushed
           down to the next lower position, the value n is put onto the
           top  of the stack, and a new repetition series is initiated.
           After the new series has been fully repeated, the value n is
           removed  from  the  stack, the value m is again put onto the
           top of the stack, and repetition of the series for m is con-
           tinued with the record following R.

       RPTO=-1
           This  value is a special value that may be used as the repe-
           tition output value with a  given  record.   If  a  previous
           record  has  been  marked as a repetition base with value 1,
           then the record in question becomes a  repetition  initiator
           with  value  1  in the usual way.  But if no previous record
           has been marked as a repetition base,  then  it  is  assumed
           that the first record of the individual record (with single-
           subject data, the first record of the data set) is a repeti-
           tion  base with value 1.  That is, before the next record is
           passed to PRED  all  preceding  records  of  the  individual
           record  (with  single-subject  data, all records of the data
           set), starting with the first such  record  and  up  to  and
           including  the  current record, will be once again passed to
           PRED (will be "repeated").  The value 1 is put onto the  top
           of the repetition stack.

           In  addition, for the repetition feature to work, it must be
           enabled at the outset.  This will be done if the repeat data
           item appears in the data set, or if RPTO is set to a nonzero
           value at ICALL=0 or ICALL=1.

  RPTON
      With each data record, RPTON may be set by PRED, and  it  conveys
      information  to NONMEM.  When the repetition output value is non-
      negative, RPTON is ignored.  Otherwise, RPTON may be  set  to  an
      integer that gives the number of times the repetition series ini-
      tiated by the data record is to be repeated.  With the  value  0,
      the  series  is  repeated once.  RPTON may be used as a left-hand
      quantity in $PRED, $PK and $ERROR blocks.

  RPTI
      RPTI is the "repetition input value".  With each data record,  it
      is  set  by NONMEM, and it conveys information to PRED.  RPTI may
      be used as a right-hand quantity in $PRED, $PK, $ERROR, $DES, and
      $AES blocks.

      RPTI=0
           The record being passed is not "being repeated".

      RPTI!=0
           The  record  being passed  to PRED is "being repeated".  The
           nonzero value is the length of  the  repetition  stack  (see
           above).

 By  default, with each pass through an individual record (with single-
 subject data, with each pass through the data set), and with any  data
 record  that  is  being  passed for the first time and is other than a
 repetition initiator, the output from PRED is used by NONMEM.  If  the
 record  is  a  repetition  initiator, NONMEM uses the output from PRED
 only when the repetition output value n has appeared on the repetition
 stack  for  the  first  time (as a result of RPTO being set to -n with
 this record) and when the record is being passed  for  the  last  time
 before  the output value is subsequently removed from the stack.  Oth-
 erwise, NONMEM ignores all output from PRED, except for the values set
 for these variables.

 For  example, in the case where a record is a repetition initiator, as
 is a subsequent record, where both records set RPTO to -1,  and  where
 both  records  set  RPTON  to 0, the first record is passed (at least)
 four times, and NONMEM uses the output from the first record  when  it
 is passed for the second time.

  PRDFL
      PRDFL  is the "PRED output control flag".  With each data record,
      it may be set by PRED, and  it  conveys  information  to  NONMEM.
      PRDFL  may  be  used  as  a  left-hand quantity in $PRED, $PK and
      $ERROR blocks.

      The PRED output control flag is used with very advanced  applica-
      tions  with the repetition feature.  With it, some of the default
      behaviour for when NONMEM pays  attention  to  PRED  output  (see
      above) - may be overridden.

       PRDFL!=1
           This  signals that the output from PRED with a passed record
           (except for the values of the repetition variables) is to be |
           ignored by NONMEM.

       PRDFL=1
           This  signals that the output from PRED with a passed record
           is to be used by NONMEM.

      With the PRED output control flag, PRED specifies when it is that
      NONMEM  is  to  use  PRED's  output.   However,  just as with the
      default behavior, where during a pass through  the  data,  NONMEM
      uses the output from a given data record once and once only, with
      each data record, PRDFL must be set to  1  once  and  once  only,
      either  when  the  record  is  passed  initially  or  when  it is
      repeated.  Moreover, as with the default behavior, PRDFL  can  be
      set  to  1 with a given record only after PRDFL has been set to 1
      with the previous data record of the individual record (for  sin-
      gle-subject data, with the previous data record of the data set).

      In addition, for the PRED output control flag to work, it must be
      enabled at the outset, i.e., PRDFL must be set to 1 at ICALL=0 or
      ICALL=1.   It may be enabled only when the repetition feature has
      also been enabled.

 Location prior to NONMEM 7: nmpr10

 (See repeti1, repeti2)                                                 |

 Help file repeti1 discusses the following files in the examples direc- |
 tory:                                                                  |

 repeat1.ctl                                                            |
 repeat1s.ctl                                                           |
 repeat1t.ctl                                                           |
 repeatf.ctl                                                            |

 REFERENCES: None.
