


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           VERBATIM CODE                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: FORTRAN statements
 CONTEXT: Abbreviated code

 USAGE:
 "      verbatim line

 DISCUSSION:
 Verbatim  code is FORTRAN code within a block of abbreviated code that
 is to be copied by NM-TRAN to the generated FORTRAN subroutine.  It is
 not  itself regarded as abbreviated code.  Such code is marked by hav-
 ing the character " as  the  first  non-blank  character.   The  "  is
 dropped  and  the  characters to the right are copied as-is to columns
 1-72 of the generated subroutine.  Lines of verbatim code can  include
 any  FORTRAN  statements: comment, declaration (e.g., INTEGER, COMMON,
 USE), I/O, CALL, assignment, continuation, DO,  GOTO,  etc.).   Lower-
 case  characters  can  be used freely, if the FORTRAN compiler permits
 them.

 NM-TRAN does make some modifications to verbatim code, as follows:

 Placement of code
      The verbatim code is adjusted as necessary to conform to  FORTRAN
      column  conventions.  Alphabetic text that follows the initial ",
      or that follows a statement number, is  moved  to  column  7,  as
      required by FORTRAN, except for:

 Comment lines
      If  the  character  that  immediately follows the initial " is !, |
      this conforms to the FORTRAN 90 syntax for  comment  lines.   The |
      line is copied unchanged.  If the character that immediately fol- |
      lows the initial " is C or c or " or *, this conforms to the FOR- |
      TRAN  77 syntax for comment lines.  The line is copied unchanged, |
      but C or c or " or * is replaced by !.  Example:
      "! this is a comment

 Continuation lines
      Fortran 77 continuation lines (non-blank in position 6)  are  not |
      permitted  with  NONMEM  7.   Instead,  the  line to be continued |
      should end with character &.  Example:                            |
      "      X=A &                                                      |
      "      +D/E                                                       |

 Replacement Rule
      In $PK, $ERROR, $DES, $AES, $INFN, and $PRED, in a line of abbre-
      viated  code  labels  of items in the data record are replaced by
      direct references to the data record  itself  (either  DATREC  or
      EVTREC,  as  appropriate).   Thus,  for  example, if such a label
      occurs as the left-hand side of  an  assignment  statement,  then
      without this rule, when the statement is executed, the value of a
      local variable having the label is  modified,  whereas  with  the
      rule, an item in the data record is modified.

 The  character  @,  immediately  before  a  label,  can  be used as an
 "escape" to prevent the label from being replaced.  If the character @
 immediately  follows the initial ", none of the labels occuring in the
 line are replaced.

 Labels occurring in common statements or as subroutine  arguments  are
 never replaced.

 Low-value Character
      Tab  characters (and other characters that are smaller than blank
      in the computer's collating sequence, such as carriage return ^M)
      are  permitted  in  verbatim  code.  With NONMEM 7, the last non- |
      blank character on the line is replaced by a space  if  it  is  a |
      low-value character.  This permits DOS-type line endings ^M.

 By  default, all verbatim code goes in the main section of code and is
 called MAIN verbatim code.  (The main section follows declarations and
 initial  executable  code  inserted by NM-TRAN.)  Within the main sec-
 tion, verbatim and abbreviated code may be freely mixed.  Each line of
 verbatim  code is positioned in the generated code after all code gen-
 erated from the preceding line of  abbreviated  code.   The  user  may
 explicitly specify a different location as follows.

 FIRST verbatim code
      Verbatim  lines  which  must  be positioned immediately after the
      declarations which are part of the normal subroutine header,  and
      prior  to  the FIRST executable statement of the subroutine, must
      precede the first line of abbreviated code and  must  start  with
      the line "FIRST.

 MAIN verbatim code
      FIRST  verbatim  code is normally terminated by the first line of
      abbreviated code.  If there is both FIRST and MAIN verbatim code,
      and/or  the main section is to start with verbatim code, the line
      "MAIN may be used to separate the FIRST and MAIN verbatim code.

 LAST verbatim code
      Verbatim lines which are to immediately precede the RETURN state-
      ment  from  the generated subroutine must follow the last line of
      abbreviated code and must be preceded by the line "LAST.

      Example:
       $ERROR
       Y=F*(1+ERR(1))
       "LAST
       " PRINT *,HH(1,1)

      This displays values after they have been  assigned,  immediately
      prior to the return.

 If  any  lines  of verbatim code are present in a block of abbreviated
 code, NM-TRAN generates USE statements  appropriate  to  the  kind  of
 abbreviated code and allows variables undefined in abbreviated code to
 be used as right-hand quantities in abbreviated code.

 Verbatim code is meant for use by expert users of NONMEM who are  able
 to understand the generated FORTRAN subroutine and check that the ver-
 batim code has the desired effect.

 REFERENCES: Guide IV, section IV.I 
