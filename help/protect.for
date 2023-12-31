


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         PROTECT FUNCTIONS                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PROTECT functions
 CONTEXT: Source code

 USAGE:
 X=..
 XLOG10=PLOG10(X)

 DISCUSSION:
 As  of  NONMEM  7.4,  a  series of routines are available that protect
 against domain violations, divide by zero, and  floating  point  over-
 flows.   Each  of  these routines start with the letter P, followed by
 the name of the mathematical operation they are to perform.  For exam-
 ple,  PLOG is the protective code routine that performs the LOG opera-
 tion.  In addition, there are first derivative (such as  PLOGD1),  and
 second  derivative (such as PLOGD2) companion routines available which
 NMTRAN uses for computing analytical derivatives. The source  code  of
 these routines are available in ..\source\PROTECT.f90.  If you wish to
 modify their behavior, then copy PROTECT.f90 to your run directory  as
 (say)  myprotect.f90, modify it, then refer to this modified code with
 $SUBROUTINES OTHER=myprotect.f90 (The OTHER option is  not  needed  if
 the P functions are coded explicitly or implicitly via $ABBR PROTECT.)

 The following protective code routines are available:

 For  all routines, if X=not a number, X is converted to machine preci-
 sion value, which is about 1.0E-15,  before performing an operation on
 it.   If  X>INFNTY (where INFNTY is approximately 1.0E+154), then X is
 converted to INFNTY before an operation is performed on it.

      PLOG(x):  returns LOG of x.  If x<SMALLZ, where SMALLZ is approx-
      imately 2.8E-103, then LOG(SMALLZ) is returned.

      PLOG10(x):   returns  LOG10  of  x.  If x<SMALLZ, where SMALLZ is
      approximately 2.8E-103, then LOG10(SMALLZ) is returned.

      PSQRT(x): returns SQRT of x.  If x=0.0d+00, then 0 is returned.

      PEXP(x): returns EXP  of  x.   If  x>40.0,  then  PEXP(100.0)  is
      returned (avoids overflow).

      PDZ(x):  returns  1/x  .   Protects  against  divide by zero.  If
      abs(x)<SMALLZ, then 1/SMALLZ is returned.

      PZR(x): returns x.  protects  against  zero.   If  abs(x)<SMALLZ,
      then SMALLZ is returned.

      PNP(x):  returns x.  Protects against non-positive.  If X<SMALLZ,
      then SMALLZ is returned.

      PHE(x): returns x. Protects against  high  exponent.   If  X>100,
      then 100 is returned. Thus PEXP(x)=EXP(PHE(x)).

      PNG(x):  returns  x.   Protects against negative.  If X<0.0, then
      0.0 is returned.

      PTAN(x): returns tan(x). Protects against returning  infinity  on
      inputs near pi/2.

      PATAN(x): returns atan(x). Protects against large intputs.

      PACOS(x),  PASIN(x):  returns acos(x), asin(x), respectively.  If
      |X| is between 1.0 and 1+10**(-08), then x is submitted as  1  or
      -1.  So, "dirty ones" are cleaned up, but values clearly beyond 1
      are allowed to trip up the function, so the user is aware of  the
      logical error in the code, and fix the issue.

 Instead of replacing various operations with protected code operations
 by hand, you can ask NMTRAN to automatically convert your code to pro-
 tected code with the following statement:

 $ABBR PROTECT

 NMTRAN will automatically replace all LOG (or DLOG) with PLOG, EXP (or
 DEXP) with PEXP, SQRT (or DSQRT) with PSQRT, / operations with *PDZ(),
 and B**E operations with PEXP(E*PLOG(B)), and so on.

 Note  that  the any P function name that is used explicitly or implic-
 itly via $ABBR PROTECT may not be used  for  a  user-defined  variable
 name.

 When  you  use $ABBR PROTECT, you will find a considerable improvement
 in estimation stability, regardless of estimation method used.

 REFERENCES: Guide Introduction_7
