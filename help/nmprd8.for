


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      RECURSIVE PRED INDICATOR                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY: I_REC=>IRECRSIV

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IRECRSIV

 DISCUSSION:

 A PRED routine may be recursive, i.e., for single-subject data, with a
 given data record, the PRED computation depends on the values of vari-
 ables  that PRED computes with any of the previous records.  PREDPP is
 such a routine.

 By default, PRED is taken to be nonrecursive.  At ICALL=0 or  1,  PRED
 can  declare  itself  to be recursive, by including the variable I_REC
 and setting it to 1.

 The setting of I_REC can be  changed  from  problem  to  problem,  but
 regardless  of the setting, with a given problem the setting is always
 taken to be 0 if the data for that problem are population data.

 If I_REC is set to 1 and the repetition feature is not used, the  SPE-
 CIAL  option  need  not  appear on the $COVARIANCE record; the special
 computation is automatically used.

 PREDPP sets I_REC to 1 at ICALL=1.

 I_REC can be set explicitly in abbreviated code.  With recursive $PRED
 abbreviated  code and with single-subject data, if I_REC is not set to
 a value in $PRED, then NM-TRAN causes I_REC to be set to 1.

 Location prior to NONMEM 7: nmprd8

 REFERENCES: None.
