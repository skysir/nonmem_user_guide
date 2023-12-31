


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $PROBLEM                              |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Identifies the start of a NONMEM Problem Specification
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $PROBLEM  [text]

 SAMPLE:
 $PROB  THEOPHYLLINE   POPULATION DATA

 DISCUSSION:
 The  $PROBLEM record identifies the start of a NONMEM problem specifi-
 cation.  The text becomes a heading for  the  NONMEM  printout.   This
 record is required.  If the problem is part of a superproblem,  $SUPER |
 must precede the $PROBLEM record.  If $SIZES is present, it must  pre- |
 cede  the  first  $SUPER or $PROBLEM record.  Otherwise, the first NM- |
 TRAN control record must be a  $PROBLEM  record.   A  $PROBLEM  record
 other than the first one marks the beginning of another problem speci-
 fication.

 The text must be contained on a single record, and only the  first  72
 characters  of text are used in the heading.  Spaces and semicolons in
 text are regarded as part of the text.  The text is optional.

 REFERENCES: Guide IV, section III.B.1 
