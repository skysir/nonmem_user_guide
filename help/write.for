


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            WRITE PRINT                             |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: FORTRAN statements
 CONTEXT: Abbreviated code

 SAMPLE:
   WRITE (50,*) 99,ID,CL,THETA(1),ETA(1)

   IF (ICALL.EQ.3) THEN
     WRITE (55,2) BIAS
   ENDIF

 DISCUSSION:

 Certain  forms of FORTRAN WRITE, PRINT, OPEN, CLOSE, and REWIND state-
 ments may be used anywhere in $PRED, $INFN, $PK, $ERROR,  $DES,  $AES,
 $AESINITIAL blocks.

   PRINT list ...
   PRINT *,list ...
   WRITE (unit,format) list ...
   WRITE (*,format) list ...
   WRITE (unit,*) list ...
   WRITE (*,*) list ...
   OPEN(unit)
   OPEN(unit,FILE=filename)
   CLOSE(unit)
   REWIND(unit)

 Single or double quotes around the filename are optional.  However, if
 the filename contains commas, semicolons, or parentheses, then it must
 be  surrounded  by  single  quotes ' or double quotes ".  Filename may
 also contain equal signs if it is enclosed in quotes.

 If the file is opened by NM-TRAN, filename may contain embedded spaces
 if  it  is  enclosed in quotes, and may contain at most 80 characters.
 If the file is opened by NONMEM, filename  may  not  contain  embedded
 spaces, and may contain at most 71 characters.

 unit 6, n or *.
        6 indicates the unit connected to the NONMEM output file
        n indicates the number of an alternative unit (40<n<2000).      |
        *  indicates  a FORTRAN system-dependent output (with most sys-
        tems, this is equivalent to using unit 6, but see  the  FORTRAN
        documentation) With OPEN, CLOSE and REWIND, the unit may not be
        6 or *.

 format
      *, 1, 2, 991, or 992.
        * indicates list-directed  output.  (FORTRAN  library  routines
        will select a format appropriate for the current run-time value
        of the quantity when displayed.)

        1 or 991 requests the built-in format specification:
            991  FORMAT (35F14.4)

        2 or 992 requests the built-in format specification:
            992  FORMAT (35E15.7)

        (Format specification numbers 991 and 992 are used in generated
        code  to  avoid conflict with statement numbers arising from DO
        WHILE statements.  The one digit versions are perhaps easier to
        remember.)

        (No  integer format specification is provided because constants
        are stored as floating point variables in generated and library
        routines.)

 list A  list  of one or more right-hand quantities, i.e., variables or
      constants that may appear on the right in abbreviated code.   May
      not  include expressions or names of abbreviated functions.  When |
      subscripts are appropriate,  these  must  be  integer  constants, |
      declared  integer  variables,  or expressions involving such con- |
      stants and variables.

      With a PRINT statement, a list may also  include  character  con-
      stants.   A  character  constant is delimited by single or double
      quotes (' or ").  As with any Fortran character constant, a  pair
      of  adjacent  delimiters  within the constant represents a single
      character.  E.g., PRINT *,'A isn''t B' will appear in the  output
      as:
      A isn't B

      In  an  initialization  or finalization block, a list may include
      elements of the OMEGA and SIGMA arrays.  In a finalization block,
      a  list  may  include  elements  of  the  standard  error  arrays
      (SETHET,SETHETR,SEOMEG and SESIGM).

      Any array whose elements may be listed individually in a WRITE or
      PRINT  statement  may also be written in its entirety, by listing
      it without any subscripts.  Specifically, one or more of the fol-
      lowing may be included in a list:

        THETA THETAFR SETHET  SETHETR
         THSIMP ETA  THSIMPR
        OMEGA(BLOCK) OMEGA(DIAG) SEOMEG(BLOCK) SEOMEG(DIAG)
        OMSIMP(BLOCK) OMSIMP(DIAG)
        SIGMA(BLOCK) SIGMA(DIAG) SESIGM(BLOCK) SESIGM(DIAG)
        SGSIMP(BLOCK) SGSIMP(DIAG)
        IIDX CNTID IIDX,CNTID

      Note that IIDX and CNTID are arrays in a module;
      (See Objective Function Value Individual).
      When  "IIDX"  or "CNTID" is listed, then only that array is writ-
      ten.  When "IIDX,CNTID" is specified, then pairs  of  values  are
      written, one pair per line, one pair for each individual record.

      The  option  BLOCK  requests  that the entire array be written in
      full symmetric form.  The option  DIAG  requests  that  only  the
      diagonal  elements  be written.  DIAG may also be coded DIAGONAL.
      Arrays of possibly different sizes (e.g., OMEGA  and  SIGMA)  may
      not  be  listed together in the same WRITE statement.  If a WRITE
      statement requests that an entire array be written, then the only
      other  items  that  may  be  listed  with the statement are other
      entire arrays.

 When entire arrays are written, the format specification with a  WRITE
 statement  must  be  *.   With  either  a WRITE or PRINT statement, an
 appropriate format is created by the NONMEM system.  With  these  for-
 mats, elements of the standard error arrays that do not exist have the
 value 1E10 and are printed as 0.1000000E+11.  (Such elements appear in
 NONMEM  output  as dots (.........).)  If the covariance step fails or
 is not requested, all elements of the standard error arrays are 0  and
 are printed as 0.0000000E+00.

 Example:  To write the omega matrix and its standard errors in symmet-
 ric form code the following.  (Only as many elements will  be  written
 as are appropriate for the dimension of omega in the problem.)

       IF (ICALL.EQ.3) WRITE (99,*) OMEGA(BLOCK),SEOMEG(BLOCK)

 A  list  may  also  include  a vector element or the entire vector (by
 listing the vector without a subscript).

 The OPEN, CLOSE, and REWIND statements are part of  the  FORTRAN  lan-
 guage.  For  details  (i.e.,  the  relationship  of  external files to
 units), see the Language Reference Manual and the Users Guide for your
 compiler.

 REFERENCES: None.
