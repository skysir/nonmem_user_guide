


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         $COVARIANCE,$COVR                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Instructions for NONMEM Covariance Step
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $COVARIANCE  [SPECIAL] [MATRIX=] [PRINT=[E][R][S]
              [COMPRESS]
              [SLOW|NOSLOW|FAST]
              [TOL=n] [ATOL=n]
              [SIGL=n] [SIGLO=n]
              [NOFCOV]
              [PARAFILE=[filename|ON|OFF] [PARAFPRINT=n]
              [THBND=n]
              [SIRSAMPLE=[numberlist] [SIRNITER=n] [SIRCENTER=n]
              [IACCEPT=x] [IACCEPTL=x]
              [SIRDF=n] [RANMETHOD=[n|S|m|P] ]
              [SIRPRINT=n] [FILE=filename] [FORMAT=s]
              [SIRTHBND=n]
              [SIRMINWT=x] [SIRMAXWT=x] [SIR_CAPCORR=x]
              [SIRSEED=n] [SIRCLOCKSEED=[0|1]]
              [SIRPARAFILE=[filename|ON|OFF] [SIRPARAFPRINT=n]
              [PRECOND=n] [PRECONDS=TOS] [PFCOND=n] [PRETYPE=n]
              [FPOSDEF=n]
              [CHOLROFF=n] [KNUTHSUMOFF=n]
              [POSDEF=n]
              [RESUME]
              [CONDITIONAL|UNCONDITIONAL]
              [OMITTED]

 SAMPLE:
 $COVARIANCE

 DISCUSSION:
 Optional.  Requests  that  the  NONMEM Covariance Step be implemented.
 This step outputs: standard errors, covariance matrix, inverse covari-
 ance  matrix,  and the correlation form of the covariance matrix.  May
 also be coded $COVR.

 If $COV is specified, then for IMP, IMPMAP, and ITS methods,  standard
 error information will be supplied for every $EST statement.  Standard
 error information for the classical methods (METHOD=0, METHOD=1)  will
 be given only if they are the last estimation method, and only if NOF-
 COV is not specified.

 OPTIONS:

 SPECIAL
      The special computation will be used in the Covariance Step  with
      a recursive PRED subroutine.  A recursive PRED subroutine is such
      that, with single-subject data, the PRED computation with a  data
      record depends on information passed to it with any of the previ-
      ous data records.  This is the default when PREDPP is used.

 MATRIX=c
      Specifies that the covariance matrix will be different  from  the
      default  (R  sup  -1 S R sup -1).  MATRIX=R requests that 2 times
      the inverse R matrix be used.  MATRIX=S requests that 4 times the
      inverse  S  matrix be used.  (R and S are two matrices from  sta-
      tistical theory,  the Hessian and Cross-Product  Gradient  matri-
      ces,  respectively.)   With  MATRIX=R the standard errors will be |
      more consistent with other nonlinear regression  software  imple- |
      mentations.  MATRIX=R should not be used with option SPECIAL.
      MATRIX  has  no  relevance  to  the Covariance step results for a
      BAYES method.

 PRINT=[E][R][S]
      Additional outputs will be  printed  besides  the  defaults.   E:
      print  the  eigenvalues  of the correlation matrix.  R: print the
      matrix .5*R.  S: print the matrix .25*S.  PRINT=R (or S)  is  not
      needed with MATRIX=R (or S).

 COMPRESS
      Covariance  Step arrays are printed in compressed format, even if
      their size is such that NONMEM would normally print them  in  the
      usual format.

 SLOW Requests  a slower method of computation.  Required when either a
      mixture model was used along with CENTERING  on  the  $ESTIMATION
      record,  or NUMERICAL was used on the $ESTIMATION record.  If not
      present, the option will be automatically supplied in  these  two
      cases.

 NOSLOW
      Requests  a  faster  method  of computation.  This is the default
      (but see SLOW).                                                   |

 FAST (NM74)                                                            |
      This is equivalent to FAST for the $EST record.  record.  If $EST |
      FAST is set, then $COV will be set to FAST, unless SLOW or NOSLOW |
      is specified.

 TOL=n (NM72)
      Tolerance.  Used only with General Nonlinear (differential  equa-
      tion solver) ADVAN's. Sets NRD=TOL.  By default, the value set on
      $SUBROUTINES record (or $TOL or TOL subroutine) is used.  If  TOL
      is  coded  on  $COVARIANCE, it overides the default.  With NONMEM
      74, this feature is deprecated.  A user-supplied  TOL  subroutine
      should test NM_STEP and set NRD accordingly.
      TOL  has  no relevance to the Covariance step results for a BAYES
      method.

 ATOL=n (NM72)
      Absolute tolerance. Used  only  with  ADVAN9,  ADVAN13,  ADVAN14,
      ADVAN15,  ADVAN16,  ADVAN17, ADVAN18, for which TOL is a relative
      tolerance. Sets ANRD=ATOL. The default is 12 (that  is,  accuracy
      is  10**(-12)).  Usually the problem runs quickly when using this
      setting.  On occasion, however, you may want to reduce ATOL (usu-
      ally  set it equal to that of TOL), and improve speeds of up to 3
      to 4 fold.

      By default, the value set on $SUBROUTINES record (or $TOL or  TOL
      subroutine)  is  used.  If ATOL is coded on $ESTIMATION, it over-
      rides the default for that step.  If ATOL is  coded  on  $COVARI-
      ANCE,  it overrides $ESTIMATION and/or the default for that step.
      With NONMEM 74, this feature is deprecated.  A user-supplied  TOL
      subroutine should test NM_STEP and set ANRD accordingly.

 SIGL=n  SIGLO=n (NM72)
      These  options  may  be  used to specify different values for the
      Covariance step.  The step size for evaluating the R matrix (cen-
      tral difference second derivative) is set to SIGL/4.  If only the
      S matrix is evaluated (central difference first derivative),  the
      step  size  is  set  to  SIGL/3.  SIGLO is the precision to which
      individual etas are optimized.  Classical NONMEM methods in  par-
      ticular  often  require  a  different significant digits level of
      evaluation (usually more stringent) during  the  Covariance  step
      than  during Estimation Step.  For example, during the Estimation
      step, NSIG=2, SIGL=6, TOL=6 may be  sufficient,  but  during  the
      Covariance  step,  you  may need SIGL=12 TOL=12 to avoid positive
      definiteness issues.
      By default, values specified on the $ESTIMATION record are used.
      (See $ESTIMATION).
      SIGL and SIGLO have no relevance to the Covariance  step  results
      for a BAYES method.

 NOFCOV (NM72)
      No  $COV  step for any classical estimation steps.  This would be
      useful if you wanted EM  estimation  analyses  performed,  and  a
      final  FOCE  analysis  performed, but did not want the program to
      spend time on standard error assessments for FOCE, which can take
      a long time relative to the other methods.

 PARAFILE=filename
      Name  of  the  "parallel file" (the parallelization profile) that
      controls parallelization (distributed computing).   Default  file
      name if not specified: parallel.pnm or parafile name specified on
      nmfe command.
      PARAFILE=ON turns on parallelization for this $COVARIANCE record.
      PARAFILE=OFF  turns  off  parallelization  for  this  $COVARIANCE
      record.

 PARAFPRINT=n (NM74)
      The print iteration intervals to the parallelization log file can
      be controlled by this option during parallelization of  the  $COV
      step.  See also $ESTIMATION record and nmfe7 command.  Default is
      PARAFPRINT=1.

 THBND=n (NM74)
      If THBND=1, any theta boundaries specified in the  $THETA  record
      causes  NONMEM to impose a non-linear transformation of the theta
      parameters so that  the  transformed  parameters  may  vary  from
      -infinity  to  infinity.  It  does this with logistic transforma-
      tions. This is suitable during the estimation step,  but  may  be
      desirable  have this off (THBND=0) for covariance assessment, and
      assess partial derivatives of the objective function with respect
      to  the thetas themselves, or some linear transformation of these
      thetas. By default THBND=1, in keeping with the behavior of  ear-
      lier  NONMEM  versions,  which  effectively  has THBND=1. Usually
      boundaries that are fairly wide will not impact how the variance-
      covariance is assessed, such as when a lower bound of 0 is given,
      but if you have very narrow boundaries set, then this can  impact
      the  assessment  of the variance-covariance of the estimates con-
      siderably, and you may wish to set THBND=0. If no lower or  upper
      bounds  are  given to thetas in $THETA record, this option has no
      impact.

 SIRSAMPLE=[numberlist] (NM74)
      Option  SIRSAMPLE  requests  the   Sampling-Importance-Resampling
      algorithm (SIR)
      (See  Guide Introduction_7, "Importance Sampling of the Variance-
      Covariance of the Parameter Estimates").
      By default SIRSAMPLE=-1, so SIR process does not occur.   SIRSAM-
      PLE  should  be set between 300 and 10000, to indicate the number
      of random samples to be generated for each of the SIRNITER itera-
      tions.  This  will  produce SIRSAMPLE importance samples, each of
      which will be placed in the raw output file.  As of nm75, SIRSAM-
      PLE  may  contain a list of integers, one for each of the SRNITER
      iterations  (use  double  quotes  around   the   list):   SIRSAM-
      PLE="1000,200,500"  If  there are fewer than SIRNITER values, the
      last value in the list will be used as the  SIRSAMPLE  value  for
      the  remaining  iterations.   Utility  programs table_quant (fre-
      quency and quantile sorting) and table_resample (resampling)  may
      be used to analyze $COV Sampling-Importance-Resampling data.

 SIRNITER=n (NM74)
      The number of times SIR sampling should occur.  Default is 1.

 SIRCENTER=n (NM74)
      Where  the sampling (proposal) density is to be centered.  On the
      first iteration, the mean of the sampling density is at the esti-
      mate.  On subsequent iterations, the mean of the sampling density
      is at the estimate (SIRCENTER=0) or at the mean  of  the  (trans-
      formed) samples of the previous iteration (SIRSAMPLE=1).  Default
      is 0.

 SIR_CAPCORR=x (NM75)
      Limit the correlation of omegas and  sigmas  generated  form  the
      proposal  density  to  have |correlation| of SIR_CAPCORR or less.
      This is similar to PSN's cap_correlation option.

 SIRMINWT=x (NM75)

 SIRMAXWT=x (NM75)
      Limit the weight range that a sample may relative to the proposal
      density,   to   avoid   extreme   values.    Default  values  are
      SIRMINWT=0.001, SIRMAXWT=1000.0

 SIRSEED=n (NM75)
      Starting seed for SIR analysis (default=11456).

 SIRCLOCKSEED=[0|1] (NM75)
      If SIRCLOCKSEED=1 (default is 0), actual starting  seed  will  be
      10000*(seconds  after  midnight)+SEED.   This  allows  a  control
      stream to produce  different  stochastic  results  for  automated
      replications,  without  the  need to modify the seed value in the
      control stream file in each replication.

 IACCEPT=x (NM74)
      The acceptance ratio acts similarly to importance sampling in  EM
      analysis.  Default is 1.

 IACCEPTL=x (NM74)
      IACCEPTL=0 (NM74) Default is 1.  The IACCEPTL option performs the
      same as listed for IACCEPTL in "Monte Carlo  Importance  Sampling
      EM".  Default is 0.

 SIRDF=n (NM74)
      The  proposal density is to be a t distribution with n degrees of
      freedom.  Default is 0, a normal density.

 RANMETHOD=[n|S|m|P] (NM74)
      See the corresponding option of the $ESTIMATION record for impor-
      tance sampling in EM analysis.

 SIRPRINT=n (NM74)
      Set  the  console print iterations interval. This does not impact
      the iterations listed in  the  raw  output  file.   Default  SIR-
      PRINT=0.

 FILE=filename (NM74)
      By  default,  the  raw  output file is whatever was listed in the
      $EST step, or root.ext, where root is the root name of  the  con-
      trol  stream  file.  You  can re-direct SIR sample listings to an
      alternative file with this option.

 FORMAT=s (NM75)
      By default, the format for the raw output file and all additional
      output  files is whatever FORMAT was specified for the $EST step.
      You can change its format with the FORMAT option.  s defines  the
      delimiter [,|s(pace)|t(ab)] followed by a Fortran format specifi-
      cation.  The default is s1PE12.5.  For more details, see the for-
      mat help item:
      (See format).
      See INTRODUCTION TO NONMEM 7, FORMAT=s1PE11.4

 SIRTHBND=n (NM74)
      As  with the deterministic covariance assessment step, by default
      the transformed parameters are sampled,  so  that  no  sample  is
      below  the  $THETA  lower bound specification, and no higher than
      the $THETA upper bound specification.  To  allow  a  boundariless
      search  in  the original theta domain, set SIRTHBND=1. You should
      also set THBND=1, so that the deterministc covariance matrix used
      as  the proposal density is also not hindered or contorted by the
      boundaries.  Default is the value of THBND, which in turn is 0 by
      default.

 SIRPARAFILE=filename
      Name  of  the  "parallel file" (the parallelization profile) that
      controls parallelization during SIR sampling (distributed comput-
      ing).   Default  file  name  if  not  specified:  parallel.pnm or
      parafile name specified on nmfe command.
      SIRPARAFILE=ON turns  on  parallelization  for  this  $COVARIANCE
      record during SIR sampling.
      SIRPARAFILE=OFF  turns  off  parallelization for this $COVARIANCE
      record during SIR sampling.

 SIRPARAFPRINT=n (NM74)
      The print iteration intervals to the parallelization log file can
      be  controlled  by  this option during parallelization of the SIR
      sampling portion of the $COV step.  See also  $ESTIMATION  record
      and nmfe7 command.  Default is SIRPARAFPRINT=1.

 PRECOND=n (NM74)
      Options PRECOND through PRETYPE affect the preconditioning of the
      R Matrix.
      (See Guide  Introduction_7,  "Preconditioning  the  R  Matrix  to
      Improve Precision and Success Rate of $COV Step").
      By  default,  PRECOND (Preconditioning cycles") is 0, and no pre-
      conditioning of the R matrix is performed.  When PRECOND=n,  then
      up  to  n  preconditioning cycles are performed.  This is used in
      combination with the PFCOND setting.

 PRECONDS=TOS (NM74)
      By default, if preconditioning is performed, it is done on Thetas
      (T), Omegas (O), and Sigmas(S). Specify PRECONDS (Preconditioning
      types)=T to do only thetas, PRECONDS=TO to  do  only  thetas  and
      omegas, etc.

 PFCON=n (NM74)
      PFCOND means 'forced' precondition cycles. Preconditioning occurs
      exactly PFCOND times, without testing if the R matrix is positive
      definite  or  not on each preconditioning cycle. On the remaining
      PRECOND-PFCOND cycles, the R matrix is tested for positive  defi-
      niteness,  and  upon  success, will terminate the preconditioning
      cycles.  Default is PFCOND=0.

 PRETYPE=n (NM74)
      By default PRETYPE(preconditioning type)=0 and the R matrix  cor-
      rector  is  V*square_root(eigenvalue). If you set PRETYPE=1, then
      the corrector is V*square_root(eigenvalue)*Vtranspose. If you set
      PRETYPE=2,  then the corrector is the correlation version of PRE-
      TYPE=1.

 FPOSDEF=n (NM74)
      By default FPOSDEF(forced positive  definite)=0.   If  FPOSDEF=1,
      then  if  the Rmatrix is not positive definite, it will be forced
      positive definite.  If PRECOND>0, this will occur after the  PRE-
      CONDth try.

 CHOLROFF=n (NM74)
      If CHOLROFF is set to 1, then one part of the R matrix evaluation
      (before inversion) will be evaluated in  the  manner  of  earlier
      versions of NONMEM.  This is strictly for comparison with earlier
      versions for diagnostic purposes.  If CHOLROFF=0, then a part  of
      the  R matrix evaluation will undergo cholesky decomposition, and
      if not positive definite, only then will  it  undergo  evaluation
      according  to earlier versions of NONMEM.  If CHOLROFF=2, then if
      cholesky decomposition fails, it will be slightly adjusted to  be
      positive  definite.  The CHOLROFF=0 is the safest option, in that
      the Cholesky decomposition (if positive definite)  provides  more
      precision  in evaluating the R matrix.  The CHOLROFF=2 setting is
      useful if you wish to increase the rate of success  in  obtaining
      an  R matrix (followed by its inverse) that could be suitable for
      SIR sampling.  If you use CHOLROFF=2, then positive  definiteness
      has been corrected before the inverse occurs, so POSDEF will have
      no additional effect.  If you use POSDEF, then CHOLROFF should be
      0 or 1.  Default is 0.

 KNUTHSUMOFF=n (NM74)
      In NONMEM 7.4, the Knuth summing method is used to allow the most
      accurate summation of individual objective function values,  even
      with large variations in values of the individual objective func-
      tion.  To turn this off, and allow a standard summation (not rec-
      ommended  except  for comparison purposes from earlier versions),
      set KNUTHSUMOFF=1. If KNUTHSUMOFF was set in the $EST  step,  but
      not  in  the  $COV  step,  the KNUTHSUMOFF value of the last $EST
      record will be used.
      Default is 0.

 POSDEF=n (NM75)
      If POSDEF=1, then if the R matrix is not  positive  definite,  it
      will  be  forced positive definite, by method of shifting the ei-
      genvalues by an amount slightly greater than  the  most  negative
      eigenvalue.  If POSDEF=2, then R matrix is made positive definite
      by making the negative valued eigenvalues approximately 1/100  of
      the least positive eigenvalue. If POSDEF=3, then R matrix is made
      positive definite by using the absolute values of the eigenvalues
      POSDEF  is  0 (no correction) by default for classical covariance
      step, and 3 for EM methods.
      The preconditioning method is a more sophisticated, but time-con-
      suming  method of making the R matrix positive definite (See PRE-
      COND).
      Default is 0.

 CONDITIONAL
      The Covariance Step is implemented, but only when the  Estimation
      Step  terminates  successfully (in this run or in a run continued
      via $MSFI).  This is the default.

 UNCONDITIONAL
      The Covariance Step is implemented regardless of how the  Estima-
      tion  Step  terminates  (in  this  run  or in a run continued via
      $MSFI).

 RESUME (NM73)
      If MSFO=msfile was specified  in  the  Estimation  Step  for  the
      FO/FOCE/Laplace  method  and  analysis was interrupted during the
      Covariance Step, then the Covariance Step may be resumed where it
      was interrupted in a subsequent problem.  Use the $MSFI record to
      specify the MSFO file of the interrupted analysis, and the RESUME
      option of the $COV record:

      $MSFI=msfile
       ...
      $COV RESUME

 OMITTED
      The Covariance Step is not implemented.

 EXAMPLE:
 $COV UNCONDITIONAL TOL=10 SIGL=10 SIGLO=11

 REFERENCES: Guide IV, section III.B.15 
 REFERENCES: Guide V, section 9.4.2 , 10.6 
