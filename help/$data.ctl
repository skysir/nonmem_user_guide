


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $DATA                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Describes the NM-TRAN data set
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $DATA  [filename|*] [(format)] [IGNORE=c1] [NULL=c2]
        [IGNORE=(list)...|ACCEPT=(list)...]
        [PRED_IGNORE_DATA]
        [NOWIDE|WIDE] [CHECKOUT]
        [RECORDS=n1|RECORDS=label]
        [LRECL=n2] [NOREWIND|REWIND]
        [NOOPEN] [LAST20=n3] [TRANSLATE=(list)]
        [BLANKOK]
        [MISDAT=r...]
        [REPL=n...]

 SAMPLE:
 $DATA       DATAFILE

 DISCUSSION:
 This  record  specifies  the data set to be used.  It is required with
 the first problem specification.  It must precede  any  other  NM-TRAN
 control  record  that refers to specific data item types.  May also be
 coded $INFILE.

 Optional with the second or  subsequent  problem  specifications.   If
 omitted,  NONMEM re-uses the data set from the previous problem (which
 will include any modifications made via transgeneration, e.g., via use
 of NONMEM's PASS, or via simulation).

 OPTIONS:

 filename
      Name  of  the  file  containing  the data set.  Must be the first
      option.  If it contains commas, semicolons, or parentheses,  then
      it  must  be  surrounded  by  single quotes ' or double quotes ".
      Filename may also contain  equal  signs  if  it  is  enclosed  in
      quotes.   If  the file is opened by NM-TRAN, filename may contain
      embedded spaces if it is enclosed in quotes, and may  contain  at
      most  80  characters.  If the file is opened by NONMEM, the file-
      name may not contain embedded spaces, and may contain at most  71
      characters.   If  filename is the same as any option of the $DATA
      record, it must be enclosed in quotes.

      * may be coded in a problem subsequent to the  first.   This  has
      the  same  effect as omitting the $DATA record (NONMEM is told to
      re-use the previous data set), but allows the CHECKOUT option  to
      be included.  With *, no other option may be included.

 (format)
      FORTRAN format specification to be used to read the data.  Format
      codes F, E, and X may be used, but not I.  When a format is  pro-
      vided,  the  label  DROP  cannot be used on the $INPUT record and
      options WIDE and NULL may not be coded.  If omitted, NM-TRAN will
      generate a suitable FORMAT specification.

 RECORDS=n1
      The number of records to be read from the NM-TRAN data set.  Com-
      ment records are not counted.   If  NM-TRAN  does  not  drop  any
      records from its data set (see IGNORE list and ACCEPT list), then
      n1 is also the number of records written to the NONMEM data  set.
      If  NM-TRAN drops records, then the total number of records writ-
      ten to the NONMEM data set is n1  minus  the  number  of  dropped
      records.
      With   NONMEM   7.5,  records  may  also  be  dropped  using  the
      PRED_IGNORE_DATA  block  of  abbreviated  code.  The  same  total
      applies to these dropped record.  See PRED_IGNORE_DATA, below.

      If  omitted,  the  records written to the NONMEM data set are all
      the records in the NM-TRAN data set up to the end-of-file (or  up
      to  a  NONMEM  FINISH  record)  minus  the  number of comment and
      dropped records.  May also be coded NRECORDS, RECS, or NRECS.

      If the option is coded as RECORDS=label, where label  is  a  data
      item  label, NM-TRAN understands the data records for the problem
      to start with the first data record of the NM-TRAN data  set  (at
      the  place  where  the file is positioned before data records are
      read; see the NOREWIND option), and to include as well, those and
      only  those  subsequent  contiguous  data records having the same
      value of the data item as does the first record.  It  counts  the
      total  number of these data records, minus any comment or dropped
      records, and puts this number in the NONMEM control file.

      In particular, the ID label may be used  (or  alternatively,  the
      option  may be coded RECORDS=IR, RECORDS=INDREC, or RECORDS=INDI-
      VIDUALRECORD).  If a label other than  ID  is  used,  the  $INPUT
      record  must  precede  the $DATA record.  If the data are single-
      subject data, the ID  data  items  used  to  determine  the  data
      records for the problem are those labeled ID (not .ID.).

      If  there  is  more  than  one problem specification with a $DATA
      record that includes an option of the  form  RECORDS=label,  then
      either  none  of  these  $DATA  records may also include a format
      specification, or all of them must include the same format speci-
      fication.  (See  records=id).

 LRECL=n2
      The number of characters in a logical record.  Needed for certain
      operating systems (e.g., IBM/CMS).

 IGNORE=c1
      [Note: The following two options, IGNORE  and  ACCEPT  allow  the
      user  to drop records from the NM-TRAN data set prior to the run.
      With NONMEM 7.5, it is also possible to  drop  records  from  the
      NONMEM  data  set  using abbreviated code, which is more flexible
      than the criteria that can be specified in the $DATA record.
      See PRED_IGNORE_DATA, below.]

      Specifies that any data record having the character c1 in  column
      1  should be ignored, i.e., these records are not included in the
      NONMEM data set.  This allows comment records to be  included  in
      the  NM-TRAN  data set.  In general, records having the character
      c1 in column 1 will be called "comment records".
      Also permitted: IGNORE='c' or IGNORE="c",  where  c  may  be  any
      character  except  space.   IGNORE=# is the default.  That is, in
      the absence of IGNORE option, any record whose first character is
      # is treated as a comment record.

      IGNORE=@  signifies  that  any  data  record having an alphabetic
      character or @ as its first non-blank character (not just in col-
      umn  1) should be ignored.  Alphabetic characters are the letters
      A-Z and a-z.  This permits a table file having header lines to be
      used as an NM-TRAN data set.

 IGNORE=(list)
      "List"  is  a  list of one or more data item labels, with logical
      operators   and    values,    of    the    form    "label=value",
      "label.EQ.value",       "label.NE.value",       "label.GT.value",
      "label.GE.value", "label.LT.value", and "label.LE.value".   (For-
      tran  90  logical operators such as '==' '/=' '<' '<=' '>' '>=' "
      may  also  be  used.)   Thus,  the   following   are   identical:
      "label=value","label==value","label.EQ.value".   With NONMEM 7.3,
      "label.NEN.value" and "label.EQN.value" are permitted.  (There is
      no  Fortran  90  operator  for  this comparison.)  If the logical
      operator is omitted, the default is "=".  With each data  record,
      the  value of the data item with the given label and the value in
      the list are compared according to the logical operator,  and  if
      result  is "true", the record is ignored, i.e. it is not included
      in the NONMEM data set (see  example  below).  Such  records  are
      called  "dropped  records".   With  "=",  "==",  "/=', ".EQ." and
      ".NE.", the value in the data record and the value  in  the  list
      are  compared as character strings. Otherwise, they are converted
      to numeric and compared numerically.   (This  is  the  case  with |
      .NEN.  and .EQN.)  This comparison is made prior to time transla-
      tion. Hence, the TIME item cannot be compared numerically  if  it
      contains non-numeric characters such as ":".

      Note:  if  the  data  file is a table file from a previous NONMEM |
      run, values that had been integers (0,1,..) in the original  data |
      file will be real values (0.0000E+00, 1.0000E+00, ...) in the ta- |
      ble file.  A comparison for equality or inequality should now  be |
      for the real value. E.g.                                          |
      IGNORE=(OCC==1.0000E+00).                                         |
      With NONMEM 7 the default format for the table file is PE11.4, as |
      in the examples above.  The FORMAT option of $TABLE may  be  used |
      to  change  this and the values used in a subsequent  IGNORE will |
      have to be changed accordingly.                                   |
      The .NEN. and .EQN.  operators  that  are  described  above  will |
      always work.

      A  data  item  label  along  with a logical operator and value is
      called a condition.  A list may contain several conditions; these
      should be separated by commas, and the list should be enclosed in
      parentheses.  Up to 100 different conditions  altogether  can  be
      specified.  IGNORE=(list) may be used with IGNORE=c, where c is a
      character.  Multiple IGNORE options with different lists  may  be
      used.   A  list may span one or more NM-TRAN records.  The use of
      "=" after IGNORE is optional, but parentheses are  required  with
      this  form  of  IGNORE.  Values may be alphabetic or numeric, and
      may optionally be surrounded by single quotes ' or double  quotes
      ".   Quotes  are  required if a value contains special characters
      such as =.  However, a value may not contain  spaces  or  commas.
      No format specification is permitted with this form of IGNORE.

      A data item type may be dropped from the NONMEM data set by means
      of the DROP or SKIP synonym on the $INPUT record,  after  records
      are  dropped  due  to  a  condition  based on the data item type.
      E.g.,
        $INPUT ... GEN=SKIP ...
        $DATA file IGNORE=(GEN='M')
      Records having GEN equal to 'M' will be dropped, and the GEN data
      item  type  will  then  be  omitted  from the NONMEM data set.  A
      dropped data item may be any alphanumeric string (without a  data
      item delimiter - a blank or a comma).

      If  there  is more than one condition, then records satisfying at
      least one of these conditions will be dropped.   In  effect,  the
      conditions  for  dropping  a  record are connected by the implied
      conjunction ".OR.".  E.g.
        IGNORE=(GEN.EQ.1,AGE.GT.60).
      Records having GEN equal to 1 or AGE greater than 60 are dropped.
      All others are accepted.

 ACCEPT=(list)
      The  ACCEPT  list  option is identical to the IGNORE list option,
      except that it specifies conditions for  acceptance  of  records.
      An ACCEPT list cannot be used together with an IGNORE list.  How-
      ever, it may be used with IGNORE=c, where c is a character.
      E.g.
        ACCEPT=(GEN.EQ.1,AGE.GT.60).
      Records having GEN  equal  to  1  or  AGE  greater  than  60  are
      accepted.  All others are dropped.

      Suppose  it  is  desired that records be dropped that satisfy the
      logical ".AND." of several conditions.  This can  be  implemented
      by  using  an  ACCEPT  list with the negations of the conditions.
      For example, suppose that records to be ignored are those  having
      GEN=1 .AND. AGE > 60.  This may be done as follows:
      ACCEPT=(GEN.NE.1,AGE.LE.60)

 PRED_IGNORE_DATA (NM75)
      Informs  NONMEM that a PRED_IGNORE_DATA pass through the data set
      is  required.    If   the   abbreviated   code   uses   variables
      PRED_IGNORE_DATA  or  PRED_IGNORE_DATA_TEST,  this option is sup-
      plied.  The  explicit  option  is  needed  if  the  use  of   the
      PRED_IGNORE_DATA  variables  occurs  in  user-written or verbatim
      code, so that NMTRAN is unaware of it.
      (See Guide Introduction_7 "Extension to $DATA  IGNORE=st  filter-
      ing")

 NULL=c2
      Specifies  that null data items in the NM-TRAN data set are to be
      replaced in the NONMEM data  set  by  the  character  c2.   E.g.,
      NULL=0 or NULL=..
      Null data items consist of a single dot (.) or consecutive commas |
      or consecutive tab characters.  c2 may be  any  character  except
      space (" ") or semicolon (";").
      Also  permitted: NULL='c' or NULL="c", where c may be any charac-
      ter.
      If this option is omitted, NM-TRAN  replaces  each  null  with  a
      space.

 NOWIDE
      Requests  that  NM-TRAN  attempt  to  limit FDATA to 80-character
      records.  Space between adjacent columns may  be  suppressed  and
      multi-line records may be generated.  This is the default.

 WIDE Requests  that  FDATA  contain  single-line  records, and that at
      least one space separate columns.  (Records in FDATA  will  never
      be  wider  than 300 characters.)  With this option, there will be
      no FINISH (FIN) record in the NONMEM data set.

 NOREWIND|REWIND
      With the first problem specification in  a  control  stream,  the
      file is positioned at its initial point so that the first NM-TRAN
      data set in the file is used.  The options  REWIND  and  NOREWIND
      apply only with a $DATA record in a subsequent problem specifica-
      tion, and only when the file named on the record is the  same  as
      that  specified for the previous problem.  When the file named on
      the record is different from  that  specified  for  the  previous
      problem,  the file is (re)positioned at its initial point so that
      the first NM-TRAN data set in the file is used.
      REWIND: Reposition the file at its intial point so that the first
      NM-TRAN data set in the file is re-used.
      (Whether  the  file  input  to NONMEM itself will be repositioned
      depends on whether this file is FDATA or  is  one  named  in  the
      $DATA  record;  see  NONMEM  Users  Guide, Part IV for a complete
      explanation.)
      NOREWIND: Leave the file at its current position so that the next
      NM-TRAN  data set in the file is used.  The $DATA record with the
      previous problem specification must  have  included  the  RECORDS
      option (or a FINISH record must have terminated the data set used
      in the previous problem), so that NM-TRAN did not read to a phys-
      ical end-of-file.  This is the default.

  CHECKOUT
      Requests  that  NONMEM implement the data checkout mode, in which
      the PRED  routine  is  not  called  and  predictions,  residuals,
      weighted residuals and the objective function are not computed.
      May also be coded CHECKDATA.  No tasks other than $TABLE or $SCAT
      can be specified.
      With NONMEM 7.5, an additional file, FDATA.csv is  produced  that
      outputs  the contents of its input data file (typically FDATA) in
      a comma delimited file format, so you can check how NONMEM inter-
      prets  the  input data.  If the REPL option of $DATA or the REPL_
      data item is used, the replicated form of the data will appear in
      FDATA.csv.

 NOOPEN
      NM-TRAN will not open the named data file.  This permits the data
      file to be created by one problem and used in a subsequent  prob-
      lem  of the same run.  May not be used with options IGNORE, DROP,
      or when data items ID, MDV, or EVID must be generated by NM-TRAN.
      With  NOOPEN,  a  format  specification is required.  No day-time
      translation takes place.

 LAST20=n3
      Override the LAST20 constant in  resource/TRGLOBAL.f90  (default: |
      50).
      One or two digit years > LAST20 are assumed to be in the 1900's,
      One or two digit years <= LAST20 are assumed to be in the 2000's.
      E.g,.  suppose LAST20=50. Then two digit years are interpreted as
      follows:
       00-50 = 2000-2050
       51-99 = 1951-1999
      LAST20=-1 can be used when two digit years span  the  year  2051.
      All  two  digit  years will be assumed to be in the same century.
      If year is recorded with four digits, it is always processed cor-
      rectly and the value of LAST20 is of no consequence.

 TRANSLATE=(list)
      "list"  describes modifications to be made to the contents of the
      data file.  It may contain one of:
      TIME/F, TIME/F/D
      and/or one of
      II/F, II/F/D

      F ("factor") may be an integer or a real value.  If F is  a  real
      number,  the  translated value in FDATA will have the same number
      of digits after the decimal point.  If F is an integer and  D  is
      omitted, there will be 2 digits.  Alternately, the number of dig-
      its may be specified explicitly by D ("digits").  If D is a  real
      number,  it  is truncated to integer.  If D is specified as 0, it
      defaults to 2.  The maximum value of D is 12.  The number of dig-
      its that may be requested in F is limited by the precision of the
      computer.

      For example, either of the following can be used to request  val-
      ues of TIME in FDATA that have 4 digits to the right of the deci-
      mal point:
      TIME/1.0000
      TIME/1/4

      Another example is
      II/0.01/6
      which divides II values by 0.01, and writes 6 digits to the right
      of the decimal point.

      If F is specified as "24" (or 24.0..), the options involving TIME
      (II) can be used to convert the units of  time  (of  the  steady-
      state  interval)  from hours to days.  The TIME (II) data item is
      first processed as if the option  were  not  present.   Then  the
      resulting value is divided by F.

      Note:  The value of TIME is divided by F, whether or not day-time
      translation occurs (i.e., whether or not relative times are being
      computed by NM-TRAN).  Similarly, the value of II is divided by F
      whether or not ":" appears in any II data value.

 BLANKOK
      Specifies that blank lines are permitted in the NM-TRAN data set.
      With all versions prior to NONMEM VI 2.0 a blank line was permit-
      ted, and was copied to the NONMEM data set. A warning message was
      issued.  With later versions, NM-TRAN stops with an error message
      when there is a blank line  in  the  NM-TRAN  data  set.   Option
      BLANKOK   restores   the   previous   behavior.    There   is  no

      abbreviation. BLANKOK must be coded in full.

 MISDAT=r (NM74)
      A numerical value indicating a missing data  value  in  the  data
      set,  which  is  displayed on $TABLE table outputs, but is safely
      interpreted as 0 by other steps of NONMEM. May be used up  to  20
      times.
      Example: $DATA mydatafile MISDAT=1.0E-99 MISDAT=1.0E-102.

 REPL=n (NM75)
      When  the REPL=n option of $DATA is coded, the NONMEM data set is
      considered to be a template data set.  NONMEM  itself  replicates
      the template data set n times at the start of the problem to cre-
      ate an expanded NONMEM data set.
      (Note that the NONMEM data set is typically the file FDATA gener-
      ated  by  NM-TRAN.  If there is nothing for NM-TRAN to change and
      the format is supplied on the $DATA record, the file named on the
      $DATA  record  is  the  NONMEM data set.  See NONMEM Users Guide,
      Part IV.)
      The REPL option may be used with the REPL_ data  item.   If  both
      are  used  the REPL_ data item applies first, and the REPL option
      applies second.
      See REPL_ data item.
      The REPL option and REPL_ data item are meant  to  be  used  with
      $SIMULATION or $DESIGN.
      For  important  information,  see  Guide  Introduction_7  Section
      "$DATA REPL".
      (See Guide Introduction_7 "Extension to $DATA  IGNORE=st  filter-
      ing")

 Another  change  in  NONMEM  VI  2.0 is that tab characters (and other
 characters that are smaller than blank  in  the  computer's  collating
 sequence,  such  as carriage return ^M) are permitted in NM-TRAN input
 files.  In the NM-TRAN data set they are treated like commas, i.e., as
 field delimiters.  In the NM-TRAN control stream they are converted to
 spaces.  They are left unchanged in verbatim code.  With NONMEM 7, the |
 last non-blank character on the line is replaced by a space if it is a |
 low-value character.

 Note: The character ":" in TIME or II  data  items  requests  day-time
 translation  of  TIME  or  II  values. These values must have the form
 hh:mm (i.e., hours:minutes).  With NONMEM 7.3, values  may  also  have |
 the form hh:mm:ss (i.e., hours:minutes:seconds).                       |
 (See Guide IV, section II.C.2)

 REFERENCES: Guide IV, section III.B.5 
 REFERENCES: Guide V, section 6.4 
 REFERENCES: Guide Introduction_7
