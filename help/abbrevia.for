


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          ABBREVIATED CODE                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: FORTRAN-like statements
 CONTEXT: Abbreviated code

 SAMPLE:
 $PK    CL=THETA(1)+ETA(1)

 DISCUSSION:
 Abbreviated  code  is  FORTRAN-like language within an NM-TRAN control
 stream that is used to specify a model.  The $PRED, $PK, $ERROR, $DES,
 $AES, $AESINITIAL, $TOL, $MIX, and $INFN records each begin a block of
 abbreviated code.  Such a block ends with  either  the  appearance  of
 another  record  beginning with a "$", or with the absence of any fur-
 ther records in the control stream.

 The general rules are:
    One statement of abbreviated per line, starting anywhere on the line.
    No statement numbers (unlike with FORTRAN).
    Comments may be included on any line after the semicolon character ";".
    No statement types other than assignment,
    IF, THEN, ELSE, ELSEIF, ENDIF, DO WHILE,
    ENDDO, CALL, WRITE, PRINT, RETURN, OPEN, CLOSE, REWIND.
    (e.g., no GOTO, READ, FORMAT.)
    With NONMEM 7.2 and higher, both lower and upper case may be used for all
    user-defined and reserved words.
    Continuation lines are permitted using the character & at the end   |
    of the line.                                                        |
    A special statement type, EXIT, is permitted.

 ASSIGNMENT STATEMENTS

   An assignment statement  (A=B)  uses  right-hand  quantities  in  an
   expression B in order to define a quantity A.

   Permissible right-hand quantities are:
   previously  defined  quantities,  including  declared variables that |
   have also been used on the left side of an assignment statement.     |
   (See $abbreviated).
   constants, e.g 1, 1.1, 3E+1 3E1, 3E01, 3E-1 3E-01, 3D+1, 3D1,  3D01,
   3D-1 3D-01
   THETA,  OMEGA  and  SIGMA elements, e.g. THETA(1), OMEGA(2,2).  When
   the two subscripts for an OMEGA or SIGMA element are the same,  only
   one subscript is needed, e.g. OMEGA(2).
   labels of data items defined in the $INPUT record
   variables  ETA(n), EPS(n), ERR(n), which are understood to have mean
   0 and variances and covariances given by elements of the  OMEGA  and
   SIGMA arrays.
   operators + - * / **
   parentheses ()
   built-in  functions:  LOG (natural log), LOG10, EXP, SQRT, SIN, COS,
   ABS, TAN, ASIN, ACOS, ATAN, INT, MIN, MAX, MOD                       |
   NONMEM functions: PHI, GAMLN                                         |
     FORTRAN functions may  have  random  variable  arguments  and  the |
     appropriate partial derivatives are computed.  The partial deriva- |
     tive of ABS(X) with respect to eta is mathematically undefined  at |
     X=0.  We are arbitrarily defining it to be dX/deta. If the predic- |
     tion depends on X, X must always be either positive  or  negative. |
     If  the  argument of GAMLN is a random variable, it must always be |
     positive.  Function PHI may have a random argument but no  partial |
     derivates are computed.                                            |

     The  INT,  MOD,  MIN,  and  MAX  functions  produce  discontinuous |
     results.  No partial derivatives are computed.  If they  are  used |
     outside  of  a simulation block and the function value affects the |
     value of the objective function, then an error in the NONMEM Esti- |
     mation Step will probably occur.  If they are needed in the compu- |
     tations of DADT (differential equations; with ADVAN 6, 8  ,9,  13, |
     14,  15,  16, 17, 18), then the funcions should be used instead in |
     the PK routine, and model event time (MTIME) parameters be used to |
     set flags for DES.                                                 |
     (See Model Time examples)                                          |
     (See Cirdadian example)
   abbreviated functions (see below)
   subscripted variables (see below)

   With previous versions of NONMEM, constants were limited to 12 char-
   acters.  With NONMEM 7 a constant may be up to 30 characters.   (The
   length is specified by constant SCO in SIZES)

   Further information about Assignment statements:

   Parentheses  may  be  nested in the right-hand side of an assignment
   statement.

   Except in certain blocks of  abbreviated  code  (i.e.  a  simulation
   block,   an   initialization   block,  a  finalization  block  or  a
   pred_ignore_data block), the variables ETA(n),  EPS(n),  ERR(n)  are
   considered to be (the base) random variables, and quantities defined
   in terms of these variables (directly or indirectly) are  themselves
   considered to be random variables.  NM-TRAN sets up computations for
   derivatives of random variables with  respect  to  the  base  random
   variables.

   Left-hand  variables  that have no reserved meaning are called user-
   defined variables.  With  previous  versions  of  NONMEM  they  were
   restricted to 8 characters in length and could not include the char-
   acter _.  With NONMEM 7 a user-defined  variable  name  consists  of
   1-20 letters (A-Z), numerals (0-9), and the character '_', beginning
   with a letter.  (The length 20 is specified by constant SD in SIZES)
   (See SIZES).

   With NONMEM 7, variables MU_1, MU_2, etc., are reserved and are used
   for mu modelling.
   See INTRODUCTION TO NONMEM 7, MU Referencing
   (See MU_Model).

   Left-hand quantities VECTRA(n), VECTRB(n), VECTRC(n) may be defined,
   and these become elements of  pre-defined  vectors  VECTRA,  VECTRB,
   VECTRC.   If  VECTRA, VECTRB, VECTRC is defined without a subscript,
   it is understood to be a simple variable.  A subscript  must  be  an |
   integer constant, or a character string that is replaced by an inte- |
   ger constant using the $ABBR REPLACE feature.  Whereas  once  it  is
   defined,  a vector element may be used as a right-hand quantity - in
   the same way as any previously defined variable -  the vector itself
   may be used only as an argument to an abbreviated function.

   (See VECTORS and ABBREVIATED FUNCTIONS, below.)

   Values may be assigned to certain special variables: RPTO, RPTON and
   PRDFL - with the repetition feature
   (See Repetition Variables);
   SKIP_ (but only in a finalization block)  -  with  superproblems  or
   subproblems
   (See SKIP).

   Variables,  constants,  and  functions  are  all of double precision
   floating point type.  (E.g., K=.5 gives K the value .5, not 0 as  it
   would be under the FORTRAN convention that K is integer type.)

   Except  in  certain  blocks  of  abbreviated code (i.e. a simulation
   block, an initialization block, or a finalization block), a label in
   an  $INPUT record may not appear on the left in an assignment state-
   ment.

   There is a reserved storage area of module NMPRD4, and when defining
   the  kth  variable  in  this  area,  it may be referenced as COM(k).
   Then, even if the variable is defined in terms of random  variables,
   it is not regarded as being a random variable; derivatives of COM(k)
   with respect to the base random variables are not computed.

 CONDITIONAL STATEMENTS

   Conditional statements have one of two forms:

   IF (condition) assignment statement

   IF (condition) THEN
      abbreviated code ...
   ELSEIF (condition) THEN
      abbreviated code ...
   ELSE
      abbreviated code ...
   ENDIF

   See below for a restriction on ELSEIF.
   Conditions may include the operators .EQ., .NE., .LE.,  .GE.,  .LT.,
   .GT.,  .AND., and .OR..  The arthmetic operators may be coded as ==,
   /=, <=, >=, <, and >, respectively.  They  may  include  expressions
   that  can  be  used  as right-hand quantities with assignment state-
   ments.  They may not include parentheses except in such expressions.
   Valid:   IF (Q.EQ.(R+C)/D) ...
   Invalid:  IF (Q.EQ.R.AND.(C.GT.D.OR.E.EQ.F)) ...

   With $PRED, $PK, and $ERROR records, a condition may test the  ICALL
   argument.   However,  if  the test is IF (ICALL.EQ.2), then ELSE may
   not be used.  Special rules apply.  For  more  detail,  see  SPECIAL
   STATEMENTS (below), and help entries for these records.

   A condition may also test certain variables defined in modules
   (See Variables_in_modules, NONMEM_modules, PREDPP_modules).

   Abbreviated code in a particular THEN or ELSE clause may not contain
   multiple definitions of the same random variable.  The code below is
   invalid:
      IF (condition) THEN
      X=...an expression involving a random variable
      X= ....
      ENDIF

   Random variables cannot be defined within nested conditionals, i.e.,
   within a conditional structure beginning with an IF  and  containing
   another  IF.   The  use  of  ELSEIF ... THEN implies a nested condi-
   tional.

   A special rule applies when random variables are defined via  condi-
   tional  statements.  If a random variable is multiply defined within
   a series of IF ... THEN structures, but all  conditions  are  false,
   then  the  value  of the random variable is set to zero.  If an ELSE
   appears, then not all conditions are false.
   Consider two cases in which the following statements  are  the  only
   ones defining TVK and K, respectively:
      IF (WT.GT.0) TVK=THETA(1)*WT
   If  the  condition is false, the non-random variable TVK retains the
   value set with the previous data record.
      IF (WT.GT.0) K=THETA(1)*WT*EXP(ETA(1))
   If the condition is false, the value of the random variable K is set
   to zero.
   NM-TRAN prints a warning message when it detects such code.

   In  $PK, $ERROR, and $PRED records, recursion code may be used in an
   explicit manner, as in this example:

      IF (WT.GT.0) THEN
        K=THETA(1)*WT*EXP(ETA(1))
      ELSE
        K=K
      ENDIF

   If the condition is false, K retains its value set with the previous
   data record.

   Recursion  code  can  be  used in $PRED, $PK, and $ERROR records for
   other purposes as well.  The following two fragments of code  illus-
   trate how one can use abbreviated code to implement recursive kinet-
   ics in $PRED.  The first example works with a single bolus dose  and
   the second example works with single or multiple bolus doses.  Simi-
   lar code can be used in $PK and $ERROR.

      K=THETA(1)*EXP(ETA(1))
      IF (TIME.EQ.0) THEN
        OLDA=AMT
        T=TIME
      ENDIF
      A=OLDA*EXP(-K*(TIME-T))
      OLDA=A
      T=TIME

      K=THETA(1)*EXP(ETA(1))
      IF (TIME.EQ.0) THEN
        A=AMT
        T=TIME
      ELSE
        A=A*EXP(-K*(TIME-T))+AMT
      ENDIF
      T=TIME

   The above forms of recursion work for recursion from one data record |
   to  the next ("inter-record" recursion).  It is also possible to use |
   recursion in a do-while loop ("intra-record", or  "do-while"  recur- |
   sion).                                                               |

   Example of a do-while recursive loop using a random variable:        |

    TERM=THETA(1)*EXP(ETA(1))                                           |
    SUM=0                                                               |
    DO WHILE(condition)                                                 |
    SUM=SUM+TERM                                                        |
    ...                                                                 |
    ENDDO                                                               |

   A product loop such as                                               |
     PROD=PROD*TERM                                                     |
   is  also  possible, as are other ways the dowhile recursive variable |
   can be used, so long as the variable appears on both  sides  of  the |
   equal sign within the DOWHILE loop: V= ... V ...                     |

 EXIT STATEMENT

   The exit statement has three forms:
     EXIT
     EXIT n
     EXIT n k

   n  is  called  the  "PRED error return code"; (also called the "PRED
   exit code.")  It must be 1 or 2.  Default is 1.
   k is the user error code.  It may be omitted; if present, it must be
   integer-valued  in the range 0-999.  With NONMEM 7.5, it may also be
   in the range 1000-9999  when  issued  during  the  Simulation  Step.
   Default is 0.

   With  all  versions  of  NONMEM, the value of k is part of the error
   message in ETEXT, which is reported in the NONMEM output report  and
   in  file  PRDERR.   With NONMEM 7.5 and later, k tells NONMEM how to
   handle EXIT statements during the Simulation Step.
   (See Simulation block).
   (See nmprd1)

   The EXIT statement causes an immediate exit from the routine and, if
   PREDPP  is being used, a subsequent immediate exit from PREDPP, with
   a return to NONMEM.  It is typically used  in  an  IF  statement  to
   avoid  further  computation  of  the  users  code when the values of
   theta/eta's set by NONMEM are inappropriate  or  would  lead  to  an
   arithmetic  exception.   If such an exit occurs during a Covariance,
   Table or Scatterplot Step, or  during  computation  of  the  initial
   value  of  the  objective  function, NONMEM will abort.  If the exit
   occurs during an Estimation  or  Initial  Estimates  Step,  NONMEM's
   action depends on the error return code value:

     n=1
     Suppose  first  that  a  search  for  eta is not underway.  Then a
     search for theta is underway, and when a search  is  underway  for
     theta,  the following happens: If the NOABORT option is used, NON-
     MEM will try to continue using different values  for  theta  (i.e.
     "theta-recovery");  otherwise,  NONMEM  aborts.  If theta-recovery
     fails, NONMEM aborts.
     Next, suppose a search for eta is underway.  NONMEM  will  try  to
     continue  using  different  values  for eta (i.e. "eta-recovery").
     Suppose eta-recovery fails.  If a search for  theta  is  not  also
     underway, NONMEM aborts.  Otherwise, the above procedure regarding
     theta-recovery applies.
     (See Guide VI, section III.K)

     n=2
     NONMEM aborts immediately.

   If NONMEM aborts, and k>0, a user message such as the  following  is
   printed in the output:
   "PK SUBROUTINE: USER ERROR CODE = k"
   This  message  is  intended  to help the user distinguish which EXIT
   statement caused NONMEM to abort when more than one  EXIT  statement
   is present in the abbreviated code.

 WRITE AND PRINT STATEMENTS

   (See write print).

 SPECIAL STATEMENTS

   Optional.   Each of these special statements is permitted in special
   blocks of abbreviated code.
   (See Initialization-Finalization block, Simulation block).
   (See Expectation block, Data_Average block).
   (See Data_ingnore block).

   CALL PASS(MODE)

   Must be coded exactly as shown.  May be used only in $PRED  initial-
   ization-finalization blocks, and in $INFN records.  If present, MODE
   becomes a reserved variable with type INTEGER, and may not  be  used
   outside the block(s).

   CALL SUPP(ie,ic)

   Must be coded as shown, with constant integer values 0 or 1 in place
   of ie and ic.  A value 1 for ie  (ic)  suppresses  output  from  the
   Estimation  (Covariance) step.  A value 0 does not suppress the out-
   put from the step.  May be used only in $PRED  initialization-final-
   ization  blocks, and in $INFN records.  The ie (ic) value remains in
   force until changed by a call to SUPP.

   CALL RANDOM(n,R)

   Must be coded as shown, with a constant integer value 1-10 in  place
   of  n.  The value is the number of the random source.  R is a random
   number from this source.  May be used only in Simulation and  Expec-
   tation  blocks.  If present, R becomes a reserved variable with type
   REAL, and may not be used outside the block(s).

   CALL SIMETA(ETA)

   Must be coded exactly as shown.  May  only  be  used  in  Simulation
   blocks.   Note  that  NM-TRAN  itself provides the minimum necessary
   call to SIMETA.  This statement is used in abbreviated code only  to
   obtain  a different value of ETA, e.g., so that the eta distribution
   may be truncated:
   DO WHILE (ETA(1).GT.5)
      CALL SIMETA(ETA)
   ENDDO

   CALL SIMEPS(EPS)

   Must be coded exactly as shown.  May  only  be  used  in  Simulation
   blocks.   Note  that  NM-TRAN provides the minimum necessary call to
   SIMEPS.  This statement is used in abbreviated code only to obtain a
   different  value of EPS, e.g., so that the distribution may be trun-
   cated (see SIMETA above).

   DO WHILE( condition )

   With NONMEM 7, DOWHILE may be used  in  all  blocks  of  abbreviated |
   code.   DOWHILE  marks the beginning of a set of statements that are
   executed repeatedly until the condition is false.  The ending of the
   set of statements is marked by the statement ENDDO.  The non-FORTRAN
   syntax DO WHILE(DATA)  is  permitted  in  $PRED  initialization  and
   finalization blocks, and in the $INFN record.
   (See Initialization-Finalization block).

   RETURN

   May  be  used in all special blocks.  RETURN statements must be used
   with caution because they by-pass certain normal  final  actions  of
   the routine.

   Another  kind  of special block is the pred_ignore_data block, which
   may be part of $PRED,$PK, or $INFN blocks.   It  uses  the  reserved
   variables  PRED_IGNORE_DATA_TEST  and  PRED_IGNORE_DATA.  By setting
   values of PRED_IGNORE_DATA=1 for certain  data  records,  it  causes
   these records to be dropped from the NONMEM data set.  Ordinary For-
   tran syntax is used in such blocks.
   (See PRED_IGNORE_DATA block)
   (See Guide Introduction_7 "Extension to $DATA IGNORE=st filtering")

 VECTORS and ABBREVIATED FUNCTIONS

   Reserved variable names VECTRA, VECTRB, VECTRC, etc. may be used for
   pre-defined  vectors.   elements.   Elements  can  be defined in any
   order.

   Vectors that are defined outside of $DES (or $AES) may not  be  used

   in  the  $DES (or $AES) code.  I.e., a vector may not be an implicit
   basic PK parameter.

   If FUNCA, FUNCB, FUNCC, etc.  is defined without a subscript, it  is
   understood  to  be a simple variable, and then it may appear as such
   as a right-hand quantity.  If FUNCx appears for the first time as  a
   right-hand  quantity, it is understood to be the reserved name of an
   abbreviated function, and then it must include  a  single  argument,
   e.g.   FUNCA(VECTRA)) or FUNCA(THETA(1)).  Complete FORTRAN code for
   an abbreviated function is supplied by  the  user  (See  abbreviated
   function).

   The  argument  of an abbreviated function may be any expression that
   can be used in an assignment statement (possibly involving  explicit
   vector elements), or it can be simply a vector.  In the latter case,
   any VECTR may be used with any FUNC, e.g.  X=FUNCA(VECTRB).

   With NONMEM VI 2.0:

   Function names may be FUNCA through FUNCC.
   Reserved argument vector names may be VECTRA through VECTRC.
   A given reserved function may have at most 9 elements in  the  argu-
   ment vector, and may appear in abbreviated code at most 9 times.

   With NONMEM 7.3:

   Function  names  may  be FUNCA through FUNCI, and constant NFUNCX in
   SIZES may be used to increase the number.
   Reserved argument vector names may be VECTRA through VECTRC.

   With NONMEM 7.4:

   Function names may be FUNCA through FUNCZ.
   Reserved argument vector names may be VECTRA through VECTRZ.
   Extended  reserved  names  are  recognized.  These  are  FUNCxy  and
   FUNCxyz, where each of x, y, z stands for an alphabetic character A-
   Z, e.g., FUNCAB or FUNCABC.  Similar  extended  reserved  names  for
   vectors are also recognized: e.g, VECTRAB or VECTRABC.
   A  reserved  function  may  appear  in  abbreviated code more than 9
   times.

   The $ABBR FUNCTION and $ABBR  VECTOR  feature  allows  the  user  to
   declare the names of user-defined functions and  the names and sizes
   of user-defined argument vectors.  Restrictions such as  the  number
   of  user functions, the number of elements in a vector, and the num-
   ber of times a given function may appear  in  abbreviated  code  are
   more flexible and may be modified by the user.
   (See $abbreviated).
   (See abbreviated function).

 PSEUDO ASSIGNMENT STATEMENTS

   Pseudo-assignment  statements  are optional and provide certain spe-
   cial instructions.  If present, they must appear (in any  order)  at
   the beginning of a record, before the rest of the essential abbrevi-
   ated code comprising that record.  They may be enclosed in parenthe-
   ses, but two or more pseudo-assignment statements enclosed by paren-
   theses must be separated by a semicolon ";".  The following  pseudo-
   statement may be used in all records except $MIX.

   COMRES=-1  :  All  quantities  defined in abbreviated code are to be
   stored locally; they are not to be stored  in  the  (global)  MODULE
   NMPRD4.
   (Variables  in  MODULE  NMPRD4 can be used for communication between
   various user-routines, and can be displayed in tables  and  scatter-
   plots.)

   See  the  entries for specific records for information about pseudo-
   assignment statements special to those records.

 RESERVED VARIABLE NAMES

   Reserved variable names include names of those variables  that  have
   special meaning for the user, as well as names of variables that are
   used internally by NM-TRAN.  Entries for $PRED, $PK,  $ERROR,  $DES,
   $AES, $AESINITIAL, and $TOL describe the reserved variable names for
   these abbreviated codes.

   A variable name that is reserved in an abbreviated code may  not  be
   used as the name of a variable defined in this code, or in any other
   abbreviated code unless the pseudo-assignment statement COMRES=-1 is
   used (in some code).  If a reserved variable name is inappropriately
   used, NM-TRAN will generate an error message.

   The following names are reserved in all abbreviated codes for inter-
   nal use by NM-TRAN:
     Names of modules and of variables in the modules
     (See Variables_in_modules).
     GETETA SIMETA SIMEPS
     COMSAV NWIND ETEXT IERPRD IERPRDU  MSEC MFIRST NETEXT
     ETA1-ETA9,ETA10-ETA99
     Annnnn Bnnnnn Cnnnnn Dnnnnn Ennnnn Fnnnnn Pnnnnn Qnnnnn
     MCnnnn MEnnnn MGnnnn MTnnnn
     (where   nnnnn   is  numeric  00000-99999,  and  nnnn  is  numeric
     0000-9999)

 SUBSCRIPTED VARIABLES

   Subscripts may be used with user-defined variables that are declared |
   to  be  arrays using the $ABBR DECLARE record, and also with certain |
   reserved variables:                                                  |

   THETA,THETAFR,SETHET,THSIMP,SETHETR,THSIMPR                          |
   OMEGA,OMEGAF,SEOMEG,OMSIMP                                           |
   SIGMA,SIGMAF,SESIGM,SGSIMP                                           |
   CNTID,IIDX                                                           |
   CDEN_                                                                |
   A, A_0 (when used in a WRITE or PRINT statement)                     |
   ETA, EPS, ERR (when used in a WRITE OR PRINT statement)              |

   Where a subscript is permitted or necessary, it may  be  an  integer |
   expression.   An  integer  expression  is  an  integer  constant,  a |
   declared integer variable that has appeared on the left  (i.e.,  has |
   been  given  a  value),  a  constant defined in SIZES, an expression |
   involving integer constants and  integer  variables  and  arithmetic |
   operators +, -, *, /, **.  It may also be a character string that is |
   replaced by an integer expression using the $ABBR  REPLACE  feature. |
   The  number of subscripts must be equal to the number of dimensions. |
   If there are two subscripts, either or  both  may  be  constants  or |
   expressions.   If  a  subscripted  variable  may be used on the left |
   side, then the subscript may also be an integer expression.  E.g.,   |

   $ABBR DECLARE INTEGER IND                                            |
    ...                                                                 |
   $PK                                                                  |
    IND=1                                                               |
    X=THETA(IND+1)                                                      |

   With reserved random variables such as
   ETA(n) EPS(n) ERR(n) A(n) A_0(n)
   integer expressions may be used only if they are listed in in  WRITE
   or  PRINT statements.  Otherwise n must be an integer constant, or a
   character string that is replaced by an integer constant  using  the
   $ABBR REPLACE feature.

 NONMEM RESERVED VARIABLES

   There  are many variables that are generally internal to NONMEM, and |
   often are not needed by users except  occasionally,  which  are  not |
   explicitly  recognized  by NMTRAN, and so cannot be used in abbrevi- |
   ated code, but must be used with verbatim code ("  at  beginning  of |
   line).   A  convenient means of accessing such variables, as well as |
   letting NMTRAN allow you to use the variables in abbreviated code is |
   to  place the MODULE definitions in an include file that begins with |
   the name NONMEM_RESERVED (case insensitive).  Insert include  state- |
   ments for such files at the beginning of the block(s) of abbreviated |
   code in which you want to use them.  The user may use any  of  these |
   variables  without  using  verbatim  code.   NMTRAN "reads" the non- |
   mem_reserved file, and remembers the variables declared in there  as |
   acceptable to use.                                                   |

   A  list  of  useful  variables  and  their  meanings  are  listed in |
    ..\guides\useful_variables.pdf.  Be careful in their  use,  as  you |
   have  the  ability to change the values of these reserved variables, |
   and this could crash the system if you change the wrong thing.       |

   NONMEM_RESERVED_GENERAL in the ..\util directory has a  few  of  the |
   quite useful variables listed.  It needs to be copied to the present |
   run directory so that NM-TRAN has access to it.  Note that there may |
   be  user-defined  variable names in the control stream that inadver- |
   tently have the same names as variables in  NONMEM_RESERVED_GENERAL, |
   which  may  cause compiler errors.  If so, make a smaller version of |
   NONMEM_RESERVED_GENERAL  with  a  name  also  starting  with   "NON- |
   MEM_RESERVED"  that contains only those reserved variables of inter- |
   est.                                                                 |

   As an example,                                                       |

   NONMEM_RESERVED_GENERAL contains                                     |

   " USE NMBAYES_INT, ONLY: ITER_REPORT,BAYES_EXTRA_REQUEST,BAYES_EXTRA |

   These variables can be used as shown in example 8:                   |

   $PK                                                                  |
   include nonmem_reserved_general                                      |
   BAYES_EXTRA_REQUEST=1                                                |
    ....                                                                |
   IF(BAYES_EXTRA==1 .AND. ITER_REPORT>=0 .AND. TIME==0.0) THEN         |
   WRITE(50,*) ITER_REPORT,ID,CL,V1,Q,V2                                |
   ENDIF                                                                |

   ITER_REPORT contains the present iteration number as reported to the |
   console  or NONMEM report file.  BAYES_EXTRA and BAYES_EXTRA_REQUEST |
   are described in example8.ctl.                                       |

 RESERVED AND OTHER FUNCTIONS

   The nonmem_reserved_general  file  contains  function  declarations,
   such  as  TFI and TFD, which are convenient functions to easily con-
   vert an integer to text ("text from integer" TFI) or  double  preci-
   sion value to text ("text from double" TFD).

   If  you wish to define your own function for use in abbreviated code
   and have the information about its proper use of arguments  be  con-
   veyed  upon  its  execution, so the compiler may detect errors, then
   one method is to package the definition of the  function  in  a  USE
   module,  such  as is done in the example User-defined Reserved Func-
   tion.

 REFERENCES: Guide IV, section IV 
 REFERENCES: Guide VI, section III.K 
 REFERENCES: Guide Introduction_7
