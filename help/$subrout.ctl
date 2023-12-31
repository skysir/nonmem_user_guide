


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            $SUBROUTINES                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Lists certain subroutine choices for the NONMEM Executable
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SUBROUTINES  [subname1 = name1] [subname2 = name2] ...
               [TOL=n] [ATOL=n] [SSTOL=n] [SSATOL=n]
               [SUBROUTINES=kind]

 SAMPLE:
 $SUBROUTINES    PRED=pred

 DISCUSSION:
 Optional.  Describes  the  choice  of  subroutines for the NONMEM exe-
 cutable (also called the NONMEM  load  module).   May  also  be  coded
 $SUBS.

 One of the following is required:
 ADVAN=name (also specifies the use of PREDPP).
 PRED=name  (specifies a user-supplied PRED routine).

 OPTIONS:

 subname=name
      Subname  is  the  entry  name of a user-supplied subroutine to be
      included in the NONMEM executable.  Name is the name  of  a  file
      containing  FORTRAN source code for the subroutine.  Name is used
      by NM-TRAN as documentation (in FREPORT)  and  for  inclusion  of
      source  code (in FSUBS).  More than one such option may be speci-
      fied.  The names need not be unique.  Name may not contain embed-
      ded spaces.
      Name  may  contain  as  many characters as fit a single line.  It
      must not start with a digit.   If  name  contains  commas,  semi-
      colons,  equal  signs,  or  parentheses, it must be surrounded by
      single quotes ' or double quotes ".

      Subname may be chosen from the following categories.

      User-supplied NONMEM routines:
           CRIT MIX PRED CRIT MIX PRIOR CONTR CCONTR USMETA SPTWO

      User-supplied PRED routine:
           PRED

      Subroutines from the PREDPP library:
            ADVAN TRANS SS

           (See ss_option).

      User-supplied PREDPP routines:
            PK ERROR MODEL DES AES TOL INFN

      Other user-supplied subroutines:
            OTHER

           OTHER specifies the name of a file that will be copied  into
           FSUBS  (e.g.,  OTHER=filename).  A subroutine or function in
           file filename might be called by a user-supplied routine  or
           by  verbatim code.  With NONMEM 7.4, it might be listed on a
           $ABBR FUNCTION record.  The OTHER option may be used  up  to
           40  times, to specify the names of up to 40 such files; each
           file may contain multiple subroutines and functions.   These
           routines must be in Fortran 90 format.

 SUBROUTINES=kind
      Specifies  the  kind  of subroutines to be included in the NONMEM
      executable ("kind" must be DP - double precision).

 TOL=n
      When PREDPP is specified with an ADVAN that requires a  TOL  rou-
      tine,  this  option  can  be  used  to  supply  a NRD ("number of
      required digits") value. "n" is an integer.  This is  a  relative
      tolerance.

      For  TOL and the options that follow (ATOL, SSTOL, SSATOL), it is
      also possible to code TOL=name to specify the name of a user-sup-
      plied TOL routine, or to include $TOL abbreviated code, either of
      which allows all these values to be assigned by  compartment.   A
      user-supplied  TOL  routine  also allows values to be assigned by
      compartment and for each NONMEM step.  See also the TOL option of
      the $COVARIANCE record.  Required.

 ATOL=n
      Specifies  the  absolute  tolerance for ADVAN9, ADVAN13, ADVAN14,
      ADVAN15, ADVAN16, ADVAN17, and  ADVAN18.   Optional.  Default  is
      1.0E-12.   See  also  the  ATOL  option  of  the  $ESTIMATION and
      $COVARIANCE records.

 SSTOL=n
      Specifies the relative tolerance  for  Steady  State  evaluation.
      Optional. Default is TOL.

 SSATOL=n
      Specifies  the  absolute  tolerance  for Steady State evaluation.
      Optional. Default is ATOL.

 REFERENCES: Guide IV, section III.B.6 
