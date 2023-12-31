


 +--------------------------------------------------------------------+
 |                                                                    |
 |                   PARTIAL DERIVATIVE INDICATORS                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY: MFIRST=>IFRSTDER,MSEC=>ISECDER,IFIRSTEM

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IFRSTDER,ISECDER,IFIRSTEM

 DISCUSSION:

  MSEC
      MSEC=1  when  NONMEM  is expecting second-partial eta-derivatives
      with the current call to PRED.
      MSEC=0 when NONMEM  is  ignoring  second-partial  eta-derivatives
      with the current call to PRED.

      Second-partial  eta-derivatives  are  never  expected  unless the
      Laplacian method is used.

  MFIRST
      MFIRST=1 when NONMEM is expecting  first-partial  eta-derivatives
      with the current call to PRED.
      MFIRST=0  when  NONMEM  is ignoring first-partial eta-derivatives
      with the current call to PRED.  This variable  is  used  only  by |
      PREDPP, not by generated FSUBS.

      Location of MFIRST and MSEC prior to NONMEM 7: rocm12

  IFIRSTEM
      With  NONMEM  7.2  and  higher, first-partial eta-derivatives are
      computed by generated code in FSUBS for classical NONMEM methods,
      but  not  for  IMP,  SAEM,  and BAYES methods.  This improves the
      speed at which the problem is evaluated.   However,  on  occasion
      such  derivatives are needed, for example, when steady state val-
      ues are to be calculated, or when stochastic  differential  equa-
      tions  are  to  be evaluated.  In such cases, insert as the first
      line in a each block of  abbreviated  code  ($PK,  $ERROR,  $DES,
      $AES, $PRED) the following line of code:
       FIRSTEM=1
      First derivatives will be evaluated for the new methods as well.

      Note  that  FIRSTEM  is  a  local  variable. In FSUBS, FIRSTEM is
      copied from IFIRSTEM in MODULE NMPRD_INT prior to being used:
        FIRSTEM=IFIRSTEM
       ...
        IF (FIRSTEM == 1) THEN
         block of first derivative code
        ENDIF

      The statement FIRSTEM=1 is inserted in generated  code  prior  to
      the  first  test of FIRSTEM.  Thus, the local variable FIRSTEM is
      changed by the user, not the global variable IFIRSTEM.

      In order to implement this feature, NM-TRAN rearranges the  order
      of  statements  in  FSUBS.  Statements that compute first-partial
      eta-derivatives are collected together into blocks of  first  de-
      rivative  code  to  be  executed  only with FIRSTEM is 1.  Option
      NOFASTDER of the $ABBREVIATED record prevents NM-TRAN from  doing
      this and restores the order of statements in FSUBS to what it was
      in previous versions.
      (See Abbreviated)

 REFERENCES: Guide VI, section III.E.4 , IV.B.2 
