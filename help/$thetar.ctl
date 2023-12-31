


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $THETAR                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives Instructions for Transforming Final Thetas
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $THETAR  Fortran statements

 SAMPLE:
 $THETAR
 THETAR=EXP(THETA)
 Or
 $THETAR
  THETAR(1:NTHETA)=EXP(THETA(1:NTHETA))
  THETAPR(1:NTHP)=EXP(THETAP(1:NTHP))

 DISCUSSION:

 The  purpose of $THETAR is to transform the final theta values for the
 NONMEM report and additional output files.  The record name  may  also
 be coded as $THR.

 In  the  above  sample  code,  it  is assumed that $THETAI was used to
 transform thetas from the natural domain to the log domain for estima-
 tion  within  NONMEM,  such  as when performing linear MU referencing.
 $THETAR causes them to be output in the natural domain as well.

 This affects the following data in the output report:
 NPARAMETR
 FINAL PARAMETER ESTIMATE for THETA.
 FIRST ORDER STANDARD ERROR OF ESTIMATE
  THETA - VECTOR OF FIXED EFFECTS PARAMETERS
 FIRST ORDER COVARIANCE MATRIX OF ESTIMATE
  THETA - VECTOR OF FIXED EFFECTS PARAMETERS

 Contents of additional output files such as .coi and .cov and .cor are
 also affected.

 The  assignment statements may be any Fortran 95 statements.  They are
 copied unchanged to subroutine SUBROUTINE  THETARSUB  in  FSUBS  (also
 found in thetair.f90).

 They  may  include  array  assigment  statements  specifying the whole
 arrays or sections of arrays.

 Arguments of the subroutine are as follows.

 THETA
      Final estimates of theta.  Input.

 THETAP
      Final estimates of thetas specified on $THETAP  records  (or,  if
      the  informative names are not used, thetas corresponding to pri-
      ors, if any).  Input.

 THETAR
      New values of THETA.  Output

 THETAPR
      New values of THETA's for priors.  Output.

 Other reserved variables that may be used are as follows:

 NTHETA
      Number of thetas to be estimated.

 NTHP Number of theta priors.

 NPROB IPROB
      These can be tested in  IF  statements  so  that  values  may  be
      assigned diffently for different problems.

 If the range is not specified, NONMEM to supply the range (which is by
 default NTHETA+NTHP).

 $THR may be used with $THI, but not necessarily.

 Another example is rescaling thetas.  E.g., suppose that
 $THETAI
 THETA=THETAI/10.

 was used to rescale thetas on input.

 Then
 $THETAR
 THETAR=THETA*10.
 can be used to report them in the original domain.

 REFERENCES: Guide Introduction_7
