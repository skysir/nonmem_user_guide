


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              INCLUDE                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING:  Causes NM-TRAN to read control stream records from a differ-
 ent file.
 CONTEXT: NM-TRAN Control Record

 USAGE:
 INCLUDE filename [ n ]

 SAMPLE:
 INCLUDE datadef

 DISCUSSION:

 This record causes NM-TRAN to read control stream records from a  dif-
 ferent  file.   The  record  must start with the characters INCLUDE or
 include.  $INCLUDE and $include are also recognized.

 OPTIONS:

 filename
      Name of the file to be read.  Required.  Filename may not contain
      embedded spaces.  If it contains commas, semicolons, or parenthe-
      ses, then it must be surrounded by  single  quotes  '  or  double
      quotes  ".   Filename  may  also  contain  equal  signs  if it is
      enclosed in quotes.

 n    Optional.  The number of copies of the file.  Default is 1.

 NMTRAN opens the named file and reads it to end-of-file.  The contents
 of  the  file  may  be any portion of an NM-TRAN control stream, e.g.,
 control records and/or abbreviated code.  After reaching  end-of-file,
 if  the  number  of copies is greater than 1, NM-TRAN rewinds the file
 and re-reads it the specified number of times.  After reaching end-of-
 file  on the final (or only) copy, NMTRAN resumes reading the original
 control stream after the include record.

 There may be more than one include  record  anywhere  in  the  control
 stream, but they may not be nested.  That is, an included file may not
 contain include records.

 If an error is detected in the NMTRAN input, the line number displayed
 is  cumulative.   The number of the line in the included file is added
 to the number of lines already read from the control  stream  and  any
 previous include files, as if the included file(s) are physically part
 of the control stream.

 EXAMPLE:

   $PROBLEM model A with data set 3
   include data3def
   include modela
   $THETA 1.3 4
   $OMEGA .04
   $SIGMA 1
   $ESTIMATION

 The file data3def contains the $INPUT and $DATA statements.  The  file
 modela  contains  the  $SUBROUTINE  statement,  $PK  block, and $ERROR
 block.  one

 Verbatim comment lines may be present in an  included  file.   If  the |
 character  that  immediately follows the initial " is !, this conforms |
 to the FORTRAN 90 syntax  for  comment  lines.   The  line  is  copied |
 unchanged.  If the character that immediately follows the initial " is |
 C or c or " or *, this conforms to the FORTRAN 77 syntax  for  comment |
 lines.  The line is copied unchanged, but C or c or " or * is replaced |
 by !.

 REFERENCES: None.
