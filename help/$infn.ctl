


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $INFN                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Marks the beginning of abbreviated code for the INFN routine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $INFN
 abbreviated code

 SAMPLE:

 $INFN
 IF (ICALL.EQ.1) CALL SUPP (1,1)

 DISCUSSION:
 The  $INFN  record is used to describe initialization processing for a
 NONMEM run, or NONMEM problem, or finalization processing for a NONMEM
 problem.  It is used with PREDPP.

 General   rules   for   abbreviated   code  are  documented  elsewhere
 (See abbreviated code). Specific rules for $INFN blocks are  described
 elsewhere.
 (See Initialization-Finalization block, Finalization example).

 RECORD ORDER:

 Follows $SUBROUTINES.

 REFERENCES: None.
