


 +--------------------------------------------------------------------+
 |                                                                    |
 |                  NEW INDIVIDUAL INDICATOR: NEWIND                  |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied PK and ERROR routines

 USAGE:
      USE PROCM_INT, ONLY: NEWIND=>PNEWIF

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: PNEWIF

 DISCUSSION:

 NEWIND
      NEWIND  is  of interest only when ICALL=2, 4, 5, or 6.  It is the
      same as the NEWIND argument passed by NONMEM to PREDPP (with  the
      first data record of the current event record).  Values are:
           NEWIND=0: First event record of the data set.
           NEWIND=1:  First  event  record of the data set, THETA value
           does not differ from value at last call  with  this  record,
           and PRED is nonrecursive (See Recursive PRED Indicator), or,
           NEWIND=1:  First  event  record  of  a subsequent individual
           record.
           NEWIND=2: Subsequent event record of an individual record.

 Location prior to NONMEM 7: procm1

 REFERENCES: Guide VI, section III.I 
