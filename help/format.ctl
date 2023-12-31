


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               FORMAT                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies the format for table files and other files
 CONTEXT: Options of NM-TRAN control file records

 USAGE:

 FORMAT  and  related  options  may  be used with the following control
 records:

 $TABLE [FORMAT=s1PE11.4(default)][LFORMAT=s][RFORMAT=s][IDFORMAT=s]
 Specifies the format for writing to the output table file.

 $ESTIMATION [FORMAT=s1PE12.5 (default)] [DELIM=s]
 Specifies the format for writing to the raw output files.

 $COVAR [FORMAT=s1PE12.5(default)] [DELIM=s]
 Specifies the format for writing to the raw output files

 $CHAIN [FORMAT=s1PE12.5 (default)] [DELIM=s]
 Specifies the format for writing to the raw output files

 $RCOV  [FORMAT=s] [DELIM=s]
 Specifies the format for reading from the input file root.cov

 $RCOVI  [FORMAT=s] [DELIM=s]
 Specifies the format for reading from the input file root.coi

 $ETAS  [FORMAT=s] [DELIM=s]
 Specifies the format for reading from the input file root.eta

 $PHIS  [FORMAT=s] [DELIM=s]
 Specifies the format for reading from the input file root.phi

 DISCUSSION:

 FORMAT=s1PE12.5

      s defines the delimiter [,|s(pace)|t(ab)] followed by a   Fortran
      format specification.  The FORMAT option is case-insensitive.

      Default for  $TABLE: FORMAT=s1PE11.4
      Default for $ESTIMATION, $COVAR, $CHAIN:
        With versions of NONMEM 7.5 and higher: FORMAT=s1PE12.5
        With versions of NONMEM before 7.5: FORMAT=s1PE11.4

      As with previous versions:

      The  first character defines the delimiter.  Delimiters are not a
      part of the Fortran language.  The delimiter is stripped from the
      Fortran  format  specification, and then used by NONMEM to modify
      the file being generated.  The delimiter may be s  for  space,  t
      for  tab,  or the comma.  A comma produces a comma delimited file
      with aligned fields (so, padded with spaces).  The delimiters may
      be upper or lower case.

      New to NONMEM 7.5:

      There are two new delimiters.  c for comma delimited file with no
      spaces, or q for comma delimited file with no spaces  and  double
      quotes   around   column   names   that   have  commas  (such  as
      "OMEGA(2,1)").
      The list of delimiters is now as follows.

      s means space
      t means tab
      , means comma
      c means comma, and spaces removed
      q means comma. Spaces removed, and quotes around entries that need them.

      Two new formats are permitted, which provide both  delimiter  and
      number format:
      FORMAT=QCSV is equivalent to FORMAT=q1PG23.16
      FORMAT=CSV is equivalent to FORMAT=c1PG23.16

      May be terminated with
      line length (e.g., FORMAT=s1PE15.8:160),
      line length with continuation marker for the end of each line (e.g., FORMAT=s1PE15.8:160c),
      and continuation marker for start of each continued line (e.g. FORMAT=s1PE15.8:160cx).
      If "c" is the letter s, it stands for "space".
      If "c" is the character &, (e.g., FORMAT=s1PE15.8:160&) and & is the last character on the line,
      it must be followed by ";" (e.g., FORMAT=s1PE15.8:160&;)

      If the user-defined format is inappropriate for a particular num-
      ber, then NONMEM uses the default format for that number.

      Most of the numeric formats are similar to those of Fortran,  but
      there are exceptions.

      The syntax for the number format is Fortran based, as follows:

      For E field: xPEw.d indicates

       w total characters to be occupied by the number (including deci-
      mal point, sign, digits, E specifier, and 2 digit magnitude),

       d digits to the right of the decimal point, and x digits to  the
      left of the decimal point.

      Examples:

      E12.5: -0.12345E+02
      2PE13.6: -12.12345E+02

      If you are outputting numbers that are less than 1.0E-99, such as
      1.2236E-102, NONMEM changes how the number  is  displayed.   With
      format 1PE12.4, Fortran would delete the "E" and display the num-
      ber as 1.2236-102.  NONMEM retains the "E" and displays one fewer
      significant  digit  to make room for the extra digit in the expo-
      nent.  The final digit is truncated (no rounding),  as follows:
      1PE12.4: 1.223E-102
      To make room for a three digit exponent, you may set  the  format
      as follows:

      xPEw.dEe

      where  e is the number of digits to be provided for the exponent.
      For example

      1PE12.4E3: 1.2236E-102
      Another example of this format:
      1PE12.4E3: -2.3456E+002

      For F field: Fw.d

      indicates

      w total characters to be occupied by the number (including  deci-
      mal point, sign and digits),

      d digits to the right of the decimal point.

      Examples:

      F10.3: -0.012
      F10.3: 234567.123

      For G field: xPGw.d

      For numbers >=0.1, will print an F field number if the value fits
      into w places showing d digits, otherwise will resort  to  xPEw.d
      format.

      For numbers <0.1, will always use xPEw.d format.

 IDFORMAT=s1PE11.4 (NM75)
      This  specifies  the format for the ID column.  By default the ID
      column has the same format  as  specified  by  FORMAT.   However,
      sometimes you wish the ID to appear as an integer, in which case,
      you may set IDFORMAT as I.  Do not include the delimiter.
      Some examples:

           IDFORMAT=I Integer value, left adjusted in the field.

           IDFORMAT=I6 Integer value, right adjusted  in  the  first  6
           characters of the field

           IDFORMAT=F6.1 Floating value, with single digit to the right
           of the decimal.

      If an improper format is given, it defaults to that of FORMAT.

 LFORMAT=s  RFORMAT=s

      An alternative format description to FORMAT is RFORMAT and  LFOR-
      MAT.

      (where R=real numbers) describes the full numeric record of a ta-
      ble, so that formats for specific columns may be specified.

      LFORMAT (where L=label) specifies the format of  the  full  label
      record of a table.

      The  formats  must be enclosed in double quotes, and (), and have
      valid Fortran format specifiers. The RFORMAT and LFORMAT  options
      can  be  repeated  if  the format specification is longer than 80
      characters. Multiple RFORMAT and LFORMAT entries will be concate-
      nated to form a single format record specification.
      For example,

      LFORMAT="(4X,A4,4(',',4X,A8))"
      RFORMAT="(F8.0,"
      RFORMAT="4(',',1PE12.5))"
      Will result in the following formats submitted to a Fortran write
      statement:
      LFORMAT=(4X,A4,4(',',4X,A8))
      for the table's label record, and RFORMAT=(F8.0,4(',',1PE12.5))
      For the table's numeric records.

      If RFORMAT and LFORMAT are given, then the FORMAT option will  be
      ignored. By default, FORMAT, RFORMAT, LFORMAT specifications will
      be passed on to the next $TABLE record in a given problem  unless
      new  ones are given. To turn off an RFORMAT/LFORMAT specification
      in a subsequent table (and therefore use FORMAT instead), set

      LFORMAT="NONE"
      RFORMAT="NONE"
      Here is an example of $TABLE statements designated in  a  control
      stream file:

      $TABLE ID TIME PRED RES WRES CPRED CWRES EPRED ERES EWRES NOAPPEND ONEHEADER
      FILE=tabstuff.TAB NOPRINT,FORMAT=,1PE15.8
      $TABLE ID CL V1 Q V2 FIRSTONLY NOAPPEND NOPRINT FILE=tabstuff.PAR
      LFORMAT="(4X,A4,4(',',4X,A8))"
      RFORMAT="(F8.0,"
      RFORMAT="4(',',1PE12.5))"

      $TABLE ID ETA1 ETA2 ETA3 ETA4 FIRSTONLY NOAPPEND NOPRINT
      FILE=tabstuff.ETA,FORMAT=";F12.4"
      LFORMAT="NONE"
      RFORMAT="NONE"

      There  is  no  NMTRAN  error  checking on the RFORMAT and LFORMAT
      records, so the user must engage in trial and error to  obtain  a
      satisfactory  table output (you should set MAXEVAL=0 or MAXEVAL=1
      for the $EST step to do a quick check, so you don't  spend  hours
      on estimation only to find the RFORMAT/LFORMAT were not appropri-
      ate).  The $MSFO option of the $ESTIMATION record should be  used
      to  produce  a Model Specification File.  Such a file can be used
      via $MSFI record.

      A word of caution. The FORMAT descriptor 1P, which means move the
      decimal point to the left by 1, will be in effect for all remain-
      ing FORMAT components. For example, in
      RFORMAT="(F8.0,37(',',1PE13.6),24(',',F7.2))"

      the F field format that follows an E field format,  in  which  1P
      was  used,  will  also have the decimal placed to the left, and a
      1.00 would appear as a 10.00. To  prevent  this  from  occurring,
      revert to no decimal shift with 0P:

      RFORMAT="(F8.0,37(',',1PE13.6),24(',',0PF7.2))"

 See INTRODUCTION TO NONMEM 7, $EST: Format of Raw Output File
 See INTRODUCTION TO NONMEM 7, IDFORMAT= I
 See INTRODUCTION TO NONMEM 7, LFORMAT, RFORMAT

 REFERENCES: Guide Introduction_7
