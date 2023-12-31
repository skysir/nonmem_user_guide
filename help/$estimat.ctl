


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         $ESTIMATION,$ESTM                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Instructions for the NONMEM Estimation Step
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $ESTIMATION
              [ABORT|NOABORT|NOHABORT]
              [ATOL=n]
              [AUTO=[0|1|2|3]]
              [BAYES_PHI_STORE=[0|1]]
              [BIONLY=[0|1]]
              [BOOTDATA=[0|1]]
              [CALPHA=n]
              [CENTERING|NOCENTERING]
              [CINTERVAL=n]
              [CITER=n | CNSAMP=n]
              [CONSTRAIN=n]
              [CTYPE=[0|1|2|3|4]]
              [DERCONT=[0|1]]
              [DF=n]
              [DFS=n]
              [EONLY=[0|1]]
              [ETABARCHECK|NOETABARCHECK]
              [ETADER=n]
              [ETASAMPLES=0 | ETASAMPLES=1]
              [ETASTYPE=0 | ETASTYPE=1]
              [FILE=filename]
              [FNLETA=n]
              [FORMAT|DELIM=s]
              [FO|NOFO]
              [FPARAFILE=[filename|ON|OFF]
              [GRD=s]
              [GRDQ=s]
              [GRID=(nr,ns,r0,r1)]
              [IACCEPT=x]
              [IACCEPTL=x]
              [IKAPPA=x]
              [INTERACTION|NOINTERACTION]
              [ISAMPEND=n]
              [ISAMPLE=n]
              [ISAMPLE_M1=n]
              [ISAMPLE_M1A=n]
              [ISAMPLE_M1B=n]
              [ISAMPLE_M2=n]
              [ISAMPLE_M3=n]
              [ISCALE_MAX=x]
              [ISCALE_MIN=x]
              [KAPPA=x]
              [KNUTHSUMOFF=n]
              [LAPLACIAN|NOLAPLACIAN]
              [LEVCENTER=[0|1]]
              [LEVWT=n]
              [LIKELIHOOD|-2LOGLIKELIHOOD]
              [LNTWOPI]
              [MADAPT=n]
              [MAPCOV=n]
              [MAPINTER=n]
              [MAPITER=n]
              [MAPITERS=[0|1]]
              [MASSRESET=n]
              [MAXEVALS=n]
              [MCETA=n]
              [METHOD=kind]
              [MSFO=filename]
              [MUM=s]
              [NBURN=n]
              [NITER=n]
              [NOCOV=[0|1]]
              [NOLABEL=[0|1]]
              [NONINFETA=[0|1]]
              [NOPRIOR=[0|1]]
              [NOSUB=[0|1]]
              [NOTITLE=[0|1]]
              [NUMDER=[0|1|2|3]]
              [NUMERICAL|NONUMERICAL]
              [NUTS_BASE=x]
              [NUTS_DELTA=x]
              [NUTS_EPARAM=n]
              [NUTS_GAMMA=x]
              [NUTS_INIT=x]
              [NUTS_MASS=[B|F|D|BD|DB|BBD|BBB]]
              [NUTS_MAXDEPTH=n]
              [NUTS_OPARAM=n]
              [NUTS_REG=x]
              [NUTS_SPARAM=x]
              [NUTS_STEPINTER=n]
              [NUTS_STEPITER=n]
              [NUTS_TERM=x]
              [NUTS_TEST=n]
              [NUTS_TRANSFORM=n]
              [OACCEPT=n]
              [OLKJDF=n]
              [OLNTWOPI]
              [OMEGABOUNDTEST|NOOMEGABOUNDTEST]
              [OMITTED]
              [OPTMAP=n]
              [ORDER=xxxf]
              [OSAMPLE_M1=n]
              [OSAMPLE_M2=n]
              [OVARF=x]
              [PACCEPT=n]
              [PARAFILE=[filename|ON|OFF]
              [PARAFPRINT=n]
              [PHITYPE=n]
              [POSTHOC|NOPOSTHOC]
              [PREDICTION]
              [PRINT=n]
              [PRIORC]
              [PSAMPLE_M1=n]
              [PSAMPLE_M2=n]
              [PSAMPLE_M3=n]
              [PSCALE_MAX=n]
              [PSCALE_MIN=n]
              [RANMETHOD=[n|S|m|P] ]
              [REPEAT|NOREPEAT]
              [REPEAT1|NOREPEAT1]
              [REPEAT2|NOREPEAT2]
              [SADDLE_HESS=n]
              [SADDLE_RESET=n]
              [SEED=n]
              [CLOCKSEED=[0|1]]
              [SELECT=n]
              [SIGDIGITS|NSIGDIGITS=n]
              [SIGL=n]
              [SIGLO=n]
              [SIGMABOUNDTEST|NOSIGMABOUNDTEST]
              [SLKJDF=x]
              [SLOW=1|SLOW=2]
              [SLOW|NOSLOW|FAST]
              [SORT|NOSORT]
              [STDOBJ=x]
              [STIELTJES]
              [SVARF=x]
              [TBLN=n]
              [THETABOUNDTEST|NOTHETABOUNDTEST]
              [THIN=n]
              [TPU=n]
              [TTDF=n]
              [ZERO=list]

 SAMPLE:
 $ESTIMATION     MAXEVAL=450  PRINT=5

 DISCUSSION:
 Optional.  Requests  that  the  NONMEM Estimation Step be implemented.
 May also be coded $ESTM or $ESTIMATE.   The  Estimation  Step  obtains
 parameter estimates.
 With  NONMEM 7, multiple Estimation Steps can be implemented in a sin-
 gle problem.  A sequence of two or more Estimation Steps  will  result
 in the sequential execution of each.  Options specified in an $ESTIMA-
 TION record will carry over to the next $ESTIMATION  record  unless  a
 new  option  is  specified.  If a particular option is not used by the
 method then the option will be ignored.  The final parameter estimates
 from an Estimation Step will be passed on as the initial estimates for
 the next Estimation Step.
 (See $ESTIMATION_record_options).

 See  also "Reserved Variables that are of Interest During the  Estima-
 tion Step", at the end of this help item.

 OPTIONS:

 ABORT
      During  the  Estimation  Step,  NONMEM  does not implement theta-
      recovery when PRED sets the error return code to  1.   (The  PRED
      error  return  code  n  is  set  by the statement "EXIT n [k]" in
      abbreviated code, or by the statement IERPRD=n  in  user-supplied
      code, or by PREDPP when it detects errors.)  This is the default.

 NOABORT
      During  the  Estimation  Step,  NONMEM implements theta-recovery,
      i.e., attempt to avoid values of theta which result in PRED error
      return  code  1.  In addition, most non-positive Hessian matrices
      will be forced to be positive definite, allowing the  program  to
      continue,  and  abnormal  termination of the estimation step will
      occur less often.

 NOHABORT
      Perform positive definite correction at all levels of the estima-
      tion.   This  can  hide  a serious ill-posed problem, so use with
      care.

 ATOL=n  (NM72)
      Absolute tolerance. Used  only  with  ADVAN9,  ADVAN13,  ADVAN14,
      ADVAN15,  ADAN16,  and ADVAN17, ADVAN18, for which TOL is a rela-
      tive tolerance. Sets ANRD=ATOL. The default is 12 (that is, accu-
      racy  is 10**(-12)).  Usually the problem runs quickly when using
      this setting.  On occasion, however, you may want to reduce  ATOL
      (usually  set  it equal to that of TOL), and improve speeds of up
      to 3 to 4 fold.

      By default, the value set on $SUBROUTINES record (or $TOL or  TOL
      subroutine)  is  used.  If ATOL is coded on $ESTIMATION, it over-
      rides the default for that step.  If ATOL is  coded  on  $COVARI-
      ANCE,  it overrides $ESTIMATION and/or the default for that step.
      With NONMEM 74, this feature is deprecated.  A user-supplied  TOL
      subroutine should test NM_STEP and set ANRD accordingly.

 AUTO=0 (NM73)
      NONMEM  does  not provide best settings of certain options.  This
      is the default.

 AUTO=1 (NM73)
      Several options will be set by NONMEM that will allow  best  set-
      tings  to  be determined.  User may still over-ride those options
      set by auto, by specifying them on the  same  $EST  record.   The
      AUTO option is ignored by the FO/FOCE/Laplace methods.

 AUTO=2 (NM74)
      AUTO=2  may be user with NUTS estimation to setup the alternative
      sampling strategy "Matt trick".

 AUTO=2 (NM74)
      AUTO=3 may be user with NUTS estimation to setup the  alternative
      sampling strategy of eta sampling.
      See INTRODUCTION TO NONMEM 7 for more information.

 BAYES_PHI_STORE=[0|1] (NM75)
      If  BAYES_PHI_STORE=1  then  phi  and  eta values from each BAYES
      iteration will be stored in root.iph.

 BIONLY=[0|1] (NM75)
      BIONLY stands for Bayesian individual parameters only,  and  when
      set  to 1, will create new samples of individual parameters only,
      but will keep the population parameters fixed.  May be used  with
      BAYES_PHI_STORE=1  or  write  statements capturing the individual
      parameters samples.  See INTRODUCTION TO NONMEM 7 for more infor-
      mation.

 BOOTDATA=[0|1] (NM75)
      By  default  (BOOTDATA=0),  when data are selected based on $SIML
      BOOSTRAP, the randomly selected subjects are analyzed during  the
      subsequent  estimation  method.  If BOOTDATA=1, then the subjects
      not selected are analyzed.
      See INTRODUCTION TO NONMEM 7 for more information.

 CALPHA=n
      n is a value between 0.01 and 0.05.  Alpha error rate to use  for
      linear   regression  test  to  assess  statistical  significance.
      Default is 0.05.

 CENTERING
      Requests that the average conditional estimates of  each  eta  be
      constrained  to  be  close to 0.  May only be used with METHOD=1.
      Not permitted with INTERACTION.

 NOCENTERING
      Requests that the average conditional estimates of each  eta  not
      be constrained.  This is the default.

 CINTERVAL=n
      Every  n  iterations is submitted to the convergence test system.
      If CINTERVAL=0, then a best CINTERVAL will be found, then used.

 CITER=n
      n is the number of latest PRINT  or CINTERVAL iterations on which
      to  perform  a linear regression test (where independent variable
      is iteration number, dependent variable is parameter value).   If
      CITER=10,  then 10 of the most recent PRINTed or CINTERVAL itera-
      tions are used for the linear regression test.   Default  is  10.
      May  also  be  coded CNSAMP.  If CINTERVAL is not specified, then
      the PRINT option is used.

 CONSTRAIN=n (NM72)
      Requests simulated annealing for parameters to slow the  rate  of
      reduction  of  the  elements of OMEGA during the burn-in phase of
      the SAEM method, allowing for a more global search of  parameters
      that minimize the objective function.  Values for n are:

      0 or 4   No simulated annealing.

      1 or 5   Requests simulated annealing for OMEGA.

      2 or 6   Requests simulated annealing for SIGMA.

      3 or 7   Requests simulated annealing for both OMEGA and SIGMA.

      Default is 1.

      When CONSTRAIN>=4, simulated annealing is also performed on diag-
      onal elements of OMEGA that are fixed  to  0  to  facilitate  any
      associated thetas.

      Simulated annealing is performed by subroutine CONSTRAINT.

      The  $ANNEAL  record facilitates EM search methods for this addi-
      tional annealing technique.  The subroutine CONSTRAINT  may  also
      be  used to provide any kind of constraint pattern on any parame-
      ters.
      (See $ANNEAL).
      The user may modify the subroutine CONSTRAINT that  performs  the
      simulated annealing algorithm.
      (See additional_output_files, raw_output_file).

 CTYPE=[0|1|2|3|4]
      CTYPE is used to define the termination test to be applied to the
      burn-in phases for SAEM and BAYES methods and to  the  estimation
      phases for the ITS, IMP, and IMPMAP methods.  (CTYPE=4 Applies to
      classical methods FO/FOCE/Laplacean.)
      CTYPE=0 indicates no termination test, the default.
      CTYPE=1 indicates that the test should be applied to  the  objec-
      tive function value, THETA's and SIGMA's but not to the OMEGA's.
      CTYPE=2  indicates  that the test should be applied to the objec-
      tive function value, THETA's, SIGMA's and  diagonal  elements  of
      OMEGA.
      CTYPE=3  indicates  that the test should be applied to the objec-
      tive function value, THETA's, SIGMA's and all OMEGA elements.
      CTYPE=4 indicates that NONMEM should test if the objective  func-
      tion  has not changed by more then NSIG digits beyond the decimal
      point over 10 iterations, even though a parameter  may  oscillate
      at  some  digit.   If this condition is satisfied, the estimation
      will terminate successfully. Applies to  FO/FOCE/Laplacean  meth-
      ods.

 DERCONT=[0|1] (NM73)
      The  derivative  continuity test (DERCONT) by default is off (0).
      When DERCONT=1, the partial derivative of the objective  function
      with  respect to thetas will perform an additional test to deter-
      mine if a backward difference assessment is more accurate than  a
      forward difference assessment.  The forward difference assessment
      can differ greatly from the  backward  difference  assessment  in
      cases  of  extreme  discontinuity  when varying certain thetas by
      even just a small amount in the model results in a  large  change
      in  objective  function,  (such  as a viral model in which a very
      small change in the potency of an  anti-viral  agent  results  in
      widely  varying  time  of return of viral load).  This results in
      standard errors being poorly assessed for thetas that do not have
      inter-subject  variances associated with them.  Setting DERCONT=1
      slows the analysis, but can provide more accurate assessments  of
      SE in such models.  The DERCONT works only for the Monte Carlo EM
      algorithms such as IMP, DIRECT, IMPMAP, and SAEM.

 DF=n The proposal density is to be a t distribution with n degrees  of
      freedom.   Default is 0, a normal density.  Used with the IMP and
      IMPMAP methods.

 DFS=n (NM73)
      Degrees of freedom for the Sigma matrix for  simulation  purposes
      by CHAIN.

      DFS=-1
           This  is  the  default.  The cholesky elements are uniformly
           varied  over  the  interval  (1-accept)*initial  value   and
           (1+accept)*initial value.

      DFS=n
           The SIGMA matrix is randomly created with an inverse Wishart
           distribution centered about the initial SIGMA  values,  with
           degrees of freedom DFS for dispersion.

      DFS=0
           As  above,  but  the  size  of  the  SIGMA matrix is used as
           degrees of freedom.

      DFS=>1000000
           SIGMA is fixed at its initial value.

 EONLY=[0|1]
      A value of 1 indicates the IMP objective function should be eval-
      uated  by  performing only the expectation step without advancing
      the population parameters. Default is 0.

 ETABARCHECK
      There is an etabar statistic (See etabar) from a  previous  prob-
      lem, and the P-value associated for the etabar statistic with the
      problem at hand relates to a hypothesis test that the true etabar
      is the same as that with the previous problem.

 NOETABARCHECK
      The  P-value  associated  for  the  etabar statistic (See etabar)
      relates to a hypothesis test that the true etabar is 0.  This  is
      the default.

 ETADER=n (NM73)
      For evaluating individual variances by numerical derivative meth-
      ods.   In  evaluating  the  MAP  objective  function,  the   term
      log(Det(V))  must  be evaluated to obtained the marginal or inte-
      grated posterior density, where V  is  the  eta  Variance  matrix
      based  on  the  subject's  posterior density. With ETADER>0, SLOW
      option may be needed.

      ETADER=0
           Expected value V, using analytical first derivatives

      ETADER=1
           Expected value V, using forward finite difference  numerical
           first  derivatives.  Needed if not all code evaluating F and
           Y derivatives with respect to eta are available for process-
           ing by NM-TRAN or in user supplied code.

      ETADER=2
           Expected  value V, using central finite difference numerical
           first derivatives.  Needed if not all code evaluating F  and
           Y derivatives with respect to eta are available for process-
           ing by NM-TRAN or in user supplied code.

      ETADER=3
           2nd derivative method of evaluating V, using numerical  sec-
           ond  derivatives  of  -log(L) with respect to etas.  This is
           equivalent to using the  "Laplace  NUMERICAL"  method,  even
           though FOCE may be selected.

 ETASAMPLES=[0|1] (NM74)
      Used with $EST METHOD=SAEM or $EST METHOD=BAYES.  ETASAMPLES=0 is

      the default.  ETASAMPLES=1 causes individual ISAMPLE  random  eta
      samples per subject, to be written to root.ets, where root is the
      root name of the control stream file.
      See "Stochastic  Approximation  Expectation  Maximization  (SAEM)
      Method" in Guide INTRODUCTION TO NONMEM 7.

 ETASTYPE=0 (NM73)
      Eta shrinkage is averaged for all subjects.  This is the default.

 ETASTYPE=1
      Eta  shrinkage  is  averaged  only among subjects that provided a
      non-zero derivative of their data likelihood with respect to that
      eta.  (See etasxi).

 FAST (NM74)
      The  FAST  option  is  available  for  FOCE/ITS methods. The FAST
      method allows use of analytical theta derivatives  to  facilitate
      FOCE  analysis.  All thetas should be MU-referenced in the manner
      described in Guide INTRODUCTION TO NONMEM  7,  "MU  Referencing".
      For  thetas  that  should  not have inter-subject variability, or
      should not be MU referenced, MU reference  it  anyway  by  adding
      addional  etas and assigning them to these thetas through MU ref-
      erencing, but set their associated omega values to 0.0 FIXED.

 FILE=filename
      Name for the raw output file.  Parameter estimates and  objective
      function  value will be printed to this file every printed itera-
      tion as indicated by the PRINT option.  Default: root.ext,  where
      root  is  the  name of the control file (not including any exten-
      sion; "nmbayes" if the name is not specified on the nmfe  command
      line).   Note  that  the names of additional output files are not
      affected by this option.  Additional output files have extensions
      .ext,  .phi, .phm, .shk, .shm .grd, .xml, .cov, .cor, .coi, .cnv,
      .smt, .rmt, .imp, .npd, .npe, .npi, .fgh, .clt, .vpd, .vpt, .ets,
      .vpt, .bfm, .iph They always have the name of the control file or
      "nmbayes" for root.
      (See additional_output_files, raw_output_file).

 FNLETA=n (NM72)
      For a thorough discussion:
      See INTRODUCTION TO NONMEM 7, General New Options for $ESTIMATION
      Record
      FNLETA=0  requests  that  the  FNLMOD  and FNLETA routines not be
      called after the Estimation and Covariance steps  are  completed.
      May  improve  run  time.   Post-hoc etas for METHOD=0 will not be
      computed.
      All table outputs (diagnostics and user selected items) will  use
      EBE's   from  final  estimation  method  (conditional  modes  for
      FO/FOCE/Laplace/ITS, conditional means for IMP/SAEM, MCMC  poste-
      rior means for BAYES).
      FNLETA=1 is the default. FNLMOD and FNLETA routines are called as
      usual.
      Diagnostics depending on EBE's such  as  CWRES,  CIWRES,  CIPRED,
      etc., will use EBE's based on the final estimation method (condi-
      tional  mode  for  FO/FOCE/Laplace/ITS,  conditional   mean   for
      IMP/SAEM,  MCMC  posterior  means for BAYES), while user selected
      items will use EBE's from the FNLETA step (eta modes).
      FNLETA=2 Requests that the estimation step not be done.  All  ta-
      ble  outputs  will  use  a  common  set of EBE's from an imported
      source.  This has value if you loaded the individual etas from an
      MSF file, or from a $PHIS/$ETAS record, and you want to calculate
      $TABLE items based on those etas, rather than from a new  estima-
      tion.
      FNLETA=3  (as  of  nm74) Like FNLETA=1, will call FNLETA, and all
      table outputs (diagnostics and  user  selected  items)  will  use
      EBE's from the FNLETA step (eta modes).

 FORMAT=s
      Format  for  the raw output file and all additional output files.
      s defines the delimiter [,|s(pace)|t(ab)] followed by  a  Fortran
      format   specification.   The  default  is  s1PE12.5.   For  more
      details, see the format help item:
      (See format).
      May also be coded DELIM.

 DELIM=s
      Same as FORMAT option.
      See INTRODUCTION TO NONMEM 7, $EST: Format of Raw Output File

 FO   Requests that the First-Order Model be  used  with  METHOD=1  and
      CENTERING.  Cannot be used with LAPLACIAN.

 NOFO Requests that the First-Order Model not be used with METHOD=1 and
      CENTERING.  This is the default.

 FPARAFILE=filename (NM74)
      Final etas (empirical Bayes estimates; EBE's) are evaluated after
      the  last  Estimation  Step (when FNLETA=1).  This computation is
      parallelized if parallelization is on for  the  final  Estimation
      Step.
      FPARAFILE=filename  specifies  a different parafile than was used
      for the Estimation Step.
      FPARAFILE=ON turns on parallelization for the EBE's.
      FPARAFILE=OFF turns off parallelization for the EBE's.
      The FPARAFILE option may be specified on any $ESTIMATION  record,
      but applies only after the last Estimation Step.

 GRD=s
      s  is  a  string  of  [G|N|D|S]'s with each symbol representing a
      THETA or SIGMA parameter in numerical order. The first m  letters
      of  GRD  refer  to the m THETA's. Then the m+1th letter refers to
      SIGMA(1,1), m+2 refers to SIGMA(2,2), etc (going along the diago-
      nal of SIGMA). Omitted symbols are assumed to be D.
      G  indicates  that the THETA should be Gibbs sampled. N indicates
      the THETA should be sampled using  the  Metropolis-Hasting  algo-
      rithm.  S indicates that the THETA is being used to model a SIGMA
      parameter.  S is used with Monte Carlo EM methods.   D  (default)
      indicates the program will decide. G and N are used only with the
      BAYES method.
      Default is DDDD...

 GRD=t1v1(n1):t2v2(n2):t3v3(n3)...
      An alternative syntax may be used.  T is parameter  type  (T  for
      theta,  S for sigma-like theta).  V is a letter (S,D or N), and n
      is a number list.  For  example,  to  specify  thetas  3,  and  5
      through  8 to be Gibbs samples, theta 4 is sigma-like, and sigmas
      1-3 are to be Metropolis-Hastings processed,
      GRD=TG(3,5-8):TS(4):SN(1-3)
      Thetas and sigmas not specified are given a  default  D  designa-
      tion.

 GRDQ=0 (NM74)
      Optional.   The gradient quick option, called GRDQ, allows thetas
      that must be gradient assessed (such as those that  are  not  mu-
      referenced) and SIGMAS to be more quickly evaluated by not evalu-
      ating the gradients for every one of the ISAMPLE random  samples,
      but chooses a subset of the most important samples.

 GRID=(nr,ns,r0,r1)
      Optional.  May  be  used with STIELTJES.  For nr, ns, r0, and r1,
      see the Introduction to Version  VI  2.0.   Briefly,  a  grid  is
      obtained  by first taking the interval [r0,r1] of the length axis
      and dividing this interval into nr equal subintervals.  ns may be
      thought  of  as  the  number  of points in a single quadrant of a
      2-dimensional  ellipse  in  n-space.   Constraints  are  nr<=100,
      ns<=9999,  0<r0<r1<1.   If r1>.9999, there is no tail region.  nr
      and  ns  should   be   integers.    The   default   values   are:
      GRID=(3,1,.6,.9).

 IACCEPT=x
      Has  different  meanings, depending on the method.  With SAEM and |
      BAYES, the scaling of OMEGA  is  adjusted  so  that  samples  are |
      accepted x fraction of the time.  See ISAMPLE_M2. Default is 0.4. |

      With Importance sampling (IMP INTERACTION), expand proposal (sam- |
      pling) density variance relative to conditional density  so  that |
      on  average conditional density/proposal density=IACCEPT (default |
      0.4).  For very sparse data or highly non-linear posterior densi- |
      ties (such as with categorical data), you may want to decrease to |
      0.1 to 0.3.                                                       |

 IACCEPT=0                                                              |
      For importance sampling only, you may set IACCEPT=0.0, and NONMEM |
      will  determine  the most appropriate IACCEPT level for each sub- |
      ject, and if necessary, will use a t- distribution  (by  altering |
      the  DF  for each subject) as well.  If IACCEPT=0, the individual |
      IACCEPT values and DF values will be listed in root.imp.          |

 IACCEPTL=x (NM74)                                                      |
      A scale to a second multi-variate normal density, to  cover  long |
      tails  in the posterior density (hence L for long tails), in com- |
      bination with the normal IACCEPT value  to  cover  the  posterior |
      density near the mode.

 IKAPPA=x] (NM74)
      Used in computing weight for individual parameters in ISAMPLE_M1B
      mode.  IKAPPA is 1.

 INTERACTION
      The dependence on etas of the model for  intra-individual  random
      error is preserved in the computation of the objective function.
      Cannot  be  used  with  CENTERING.   With NONMEM 7.3, This is the |
      default with EM/Bayes methods and is supplied if NOINTERACTION is |
      not  specified  by the user.  With NONMEM 7.4, INTERACTION is not |
      supplied if LIKELIHOOD is present.

 NOINTERACTION
      Always set etas to 0 during the  computation  of  the  model  for
      intraindividual  random  error.   This  is  the default with non-
      EM/Bayes methods.

 ISAMPEND=n (NM73)
      For SAEM, if ISAMPEND is specified  as  an  upper  integer  value
      (usually  10),  then  NONMEM will perform a ISAMPLE preprocess to
      determine the best ISAMPLE value. See also STDOBJ.

 ISAMPLE=n
      When used with the IMP or IMPMAP methods n is the number of  ran-
      dom  samples  per subject used for the expectation step.  Default
      is 300.  When used with the SAEM or BAYES method n is the  number
      of  chains used by the Metropolis-Hastings algorithm for individ-
      ual parameter estimation. The default is 2 for  SAEM  and  1  for
      BAYES.

      A  kernel  is  the  Metropolis-Hastings sampling and testing mode
      used.  The ISAMPLE_Mx options define how many times  to  generate
      and test a sample for goodness-of-fit using a given kernel. ISAM-
      PLE does not refer to a kernel, but defines the number of  chains
      that  are  maintained, each chain having their own sample genera-
      tion and testing sequence using the various kernels.  Each  chain
      retains a final sample for each subject, at the end of each iter-
      ation.

 ISAMPLE_M1=n
      n is the number of mode 1 iterations for  the  Metropolis-Hasting
      algorithm  for estimating individual parameters using the popula-
      tion means and variances as proposal density. Used with the  SAEM
      and BAYES methods.  Default is 2.

 ISAMPLE_M1A=n(NM72)
      n  is the number of mode 1A iterations for the Metropolis-Hasting
      algorithm, testing model parameters from other subjects as possi-
      ble values.  Used with the SAEM and BAYES methods.  Default is 0.

 ISAMPLE_M1B=n (NM74)
      n  is the number of mode 1B iterations for the Metropolis-Hasting
      algorithm Default is 2.

 ISAMPLE_M2=n
      n is the number of mode 2 iterations for  the  Metropolis-Hasting
      algorithm  for estimating individual parameters using the current
      parameter vector position as mean and a scaled variance of  OMEGA
      as variance. Used with the SAEM and BAYES methods.  Default is 2.

 ISAMPLE_M3=n
      n  is  the number of mode 3 iterations for the Metropolis-Hasting
      algorithm for estimating individual parameters in  which  samples
      are  generated  for each parameter separately. Used with the SAEM
      and BAYES methods.  Default is 2.

 ISCALE_MIN=x ISCALE_MAX=x (NM72)
      In importance sampling, the scale factor used to vary the size of
      the variance of the proposal density in order to meet the IACCEPT
      condition is  by  default  bounded  by  ISCALE_MIN  of  0.1,  and
      ISCALE_MAX=10.0.   Defaults  are  (1.0E-06,1.0E+06) for MCMC sam- |
      pling.  These options allow the scale factor boundary to be modi-
      fied.

 KAPPA=x (NM74)
      Used with NUTS method.  Default is 1.

 KNUTHSUMOFF=n] (NM74)
      The  Knuth summing method is used to allow the most accurate sum-
      mation of individual objective function values, even  with  large
      variations  in  values  of the individual objective function.  To
      turn this off, and allow a standard  summation  (not  recommended
      except for comparison purposes from earlier versions), set KNUTH-
      SUMOFF=1. With KNUTHSUM algorithm on by default, the SORT  option
      is  not  necessary.  Default is 0.  May also be set with $COVARI-
      ANCE record.

 LAPLACIAN
      Use the  Laplacian  method,  in  which  second  derivatives  with
      respect  to  eta  are  used.   Laplacian  may  not  be  used with
      METHOD=0.  It may be used with the  EM/Monte  Carlo  methods,  in
      which  case  the Laplacian option will be properly utilized, such
      as during MAP estimation used during IMP,  IMPMAP,  and  ITS,  or
      ignored, such as during SAEM or BAYES.
      Cannot  be  used  with  $ABBREVIATED  DERIV2=NO  unless NUMERICAL
      option is also specified.

 NOLAPLACIAN
      Do not use the Laplacian method.  This is the default.

 LEVCENTER=[0|1] (NM75)
      There is no default.  Required with $LEVEL and  $ESTIMATION.   If
      LEVCENTER=1,  this ensures the etas of super ID random levels sum
      to 0. In earlier versions of NONMEM, this was  the  default  (and
      only)  action.  To  obtain similar results as earlier versions of
      NONMEM, set LEVCENTER=1.
      If LEVCENTER=0, level etas are not forced to sum to 0.
      See INTRODUCTION TO NONMEM 7 for more information.

 LEVWT=n (NM74)
      This option applies when $LEVEL record is present.   By  default,
      LEVWT=0, and weights each level value equally, regardless of num-
      ber of subjects per level value. To weight according to number of
      subjects for that value, set LEVWT=1.

 LIKELIHOOD
      This  is  designed mainly, but not exclusively, for use with non-
      continuous observed responses ("odd-type data").  Indicates  that
      Y (with NM-TRAN abbreviated code) or F (with a user-supplied PRED
      or ERROR code) will be set to a (conditional)  likelihood.   Upon
      simulation  it  will be ignored, and the DV data item will be set
      directly to the simulated value  in  abbreviated  or  user  code.
      Also etas, if any, are understood to be population etas.  Epsilon
      variables and the $SIGMA record may not be  used.   The  L2  data
      item  may not be used.  The CONTR and CCONTR options of the $SUB-
      ROUTINES record may not be used.  NONMEM cannot obtain  the  ini-
      tial  estimate  for omega.  If the data are population, and MAXE- |
      VALS=0 is not coded, then  NOINTERACTION  is  required.   Compare
      with PREDICTION option.

 -2LOGLIKELIHOOD
      Indicates  that  Y  (with  NM-TRAN abbreviated code) or F (with a
      user-supplied PRED or ERROR code) is a -2 log (conditional) like-
      lihood.   All  remarks  for  LIKELIHOOD apply.  May also be coded
      -2LLIKELIHOOD.  Compare with PREDICTION option.

 MADAPT=n (NM74)
      Used with NUTS method.  Default is -1.

 MAPCOV=1 (NM74)
      Unused. The default is 1.

 MAPINTER=n (NM72)
      Every nth iteration, the MAP estimation should be used to provide
      parameters to the sampling density.  Thus, if MAPITER=20 and MAP-
      INTER=5, then for the first  20  iterations,  MAP  estimation  is
      used,  and  thereafter, every 5th iteration the MAP estimation is
      used.  If MAPINTER=-1, then mapinter will be turned  on  only  if
      the objective function increases consistently over several itera-
      tions.  Setting any of the above parameters to  -100  will  force
      NONMEM to select the default value for that parameter.

 MAPITER=n (NM72)
      The first n iterations are to use MAP estimation to assess param-
      eters for the sampling density.  After these  n  iterations,  the
      conditional  means  and  variances  of the pervious iteration are
      used for the sampling density parameters of  the  present  itera-
      tion.  If MAPITER=0, then the first iteration will rely on condi-
      tional means and variances that are in memory.   These  may  have
      come from an MSF file, or from a previous estimation step.

 MAPITERS=[0|1] (NM75)
      By  default,  no  MAP estimation is performed with SAEM or BAYES.
      To get good individual parameter values near the mode of the pos-
      terior  density  for  the  first  iteration  of SAEM, you can set
      MAPITERS=1.  Alternatively, you can insert the record:
      $EST METHOD=ITS NITER=0
      Followed by
      $EST METHOD=SAEM
      or
      $EST METHOD=BAYES

 MASSRESET=n (NM74)
      Mass matrix information accumulation for NUTS method.  Default is
      -1.

 MAXEVALS=n
      Maximum allowable number of evaluations of the objective function
      during the Estimation Step.  Default: a generous  number.   (Each
      evaluation  of  the  objective function requires one pass through
      the data set.  This is also referred to as  a  "function  evalua-
      tion.")   MAXEVALS=-1  may  be  specified  when a $MSFI record is
      present.  It requests that NONMEM re-use the value from the  pre-
      vious run, and is the default with $MSFI.

      MAXEVALS=0 requests that the Estimation Step be omitted.  This is
      useful, for example, with POSTHOC (see above).

 MCETA=n (NM73)

      MCETA=0
           Eta=0 is initial setting for MAP estimation  (eta  optimiza-
           tion).  This is the default.

      MCETA=1
           ETA=values  of previous iteration is initial setting for MAP
           estimation, or ETA=0, whichever gives lower objective  func-
           tion.

      MCETA>1
           MCETA-1 Random samples of ETA, using normal random distribu-
           tion with variance OMEGA, are tested.  Plus previous ETA  is
           tested,  and  ETA=0  is  tested.  Whichever gives the lowest
           objective function is used as initial setting  for  the  MAP
           optimization.

 METHOD=kind
      Values for kind are:

      0 or ZERO
           Always set etas to 0 during the computation of the objective
           function.  Also called the "first order (FO) method."   This
           is the default.

      1 or CONDITIONAL
           Use  conditional  estimates for the etas during the computa-
           tion of  the  objective  function.   METHOD=1  (without  the
           LAPLACIAN  option)  is  also  called the "first order condi-
           tional estimation (FOCE) method."  The conditional estimates
           of  the  etas are referred to as Conditional Parametric Etas
           (CPE).

           METH=COND NOLAPLACIAN is referred to as the FOCE method.
           METH=COND LAPLACE is referred to as the Laplace method.
           METH=COND NOLAPLACE CENTERING is referred to as the Centering FOCE method.
           METH=COND LAPLACE CENTERING is referred to as the Centering Laplace method.

      HYBRID
           Use conditional estimates for the etas during  the  computa-
           tion  of the objective function, with the exception of those
           etas listed in the ZERO option.  Cannot be used with  LAPLA-
           CIAN or CENTERING.

      The  following  methods  are  new to NONMEM 7.  When any of these
      methods are used, the data are inferred  to  be  population,  and
      METHOD=1  is  supplied  if  it is not already present.  The first
      four methods are referred  to  as  EM  (Expectation-Maximization)
      Methods.

      ITS  Use  the  iterative two-stage method.  This method evaluates
           the conditional mode and first order  approximation  of  the
           conditional  variance  of parameters of individuals by maxi-
           mizing the posterior density.  This integration step is  the
           same  as used in the FOCE method.  Population parameters are

           updated from individuals' conditional  mode  parameters  and
           their approximate variances by single iteration maximization
           steps.

      IMP  Use the Monte-Carlo Importance  Sampling  Expectation  Maxi-
           mization method.  This method evaluates the conditional mean
           and variance of parameters of  individuals  by  Monte  Carlo
           sampling.  It uses the posterior density, which incorporates
           the likelihood of parameters relative  to  population  means
           and variances with the individual's observed data.  The nor-
           mal density near the mean or mode of the posterior  is  used
           as a proposal density, then weighted according to the poste-
           rior density as a correction.

      IMPMAP
           Use the Importance Sampling method assisted by Mode a Poste-
           riori  (MAP)  estimation.   At  each  iteration, conditional
           modes and conditional first order variances are evaluated as
           in  the  ITS or FOCE method.  These are then used as parame-
           ters to the multivariate normal  proposal  density  for  the
           Monte-Carlo importance sampling step.

      SAEM Use  the  Stochastic  Approximation Expectation Maximization
           method.  As in importance sampling, random samples are  gen-
           erated  from normal proposal densities.  However, instead of
           always being centered at the mean or mode of  the  posterior
           density,  the  proposal  density is centered at the previous
           sample position.

      BAYES
           Use the Markov Chain Monte Carlo  (MCMC)  Bayesian  Analysis
           method.  The goal of the MCMC Bayesian analysis is to obtain
           a large sample set of probable population parameters.  Vari-
           ous summary statistics of the population parameters may then
           be obtained such as means and confidence ranges.

      DIRECT
           Requests Monte Carlo Direct  Sampling.   Creates  completely
           independent samples (unlike MCMC), and there is no chance of
           causing bias if the sampling density is not  similar  enough
           to  the  conditional  density  (unlike IMP).  However, it is
           very inefficient, requiring ISAMPLE  settings  of  10000  to
           300000 to properly estimate the problem.

      NUTS (NM74)
           Requests  No U-Turn Sampling (NUTS) Markov Chain Monte Carlo
           (MCMC) Bayesian Analysis Method.   Options  unique  to  this
           method  are  listed  alphabetically  under  NUTS_....  Other
           options of interest with their defaults are as follows:
                MASSRESET=-1
                MADAPT=-1
                KAPPA=1
                TTDF=0
                OLKJDF=0
                OVARF=1
                SLKJDF=0
                SVARF=1

      CHAIN
           Allows the user to create a series of random initial  values
           of  THETAs and OMEGA's, or for reading in initial population
           parameters from a file of rectangular (rows/column)  format.
           Applies only to the Estimation Step.

 LNTWOPI (NM74)
      The  objective function is reported including the N*LOG(2pi) con-
      stant term, where N is the total number of  normally  distributed
      data  values  in  the data set.  Compare OLNTWOPI. Either or both
      may be used.

 MSFO=filename
      A Model Specification File is output to a  file  with  the  given
      filename.  Filename may not contain embedded spaces.  If filename
      contains commas, semicolons, or parentheses, then it must be sur-
      rounded  by  quotes  ('  or  ").  Filename may also contain equal
      signs if it is enclosed in quotes.  Filename may contain at  most
      71  characters.   If  filename  is  the same as any option of the
      $ESTIMATION record, it must be enclosed in quotes.  If the  $NON-
      PARAMETRIC  record is present and also specifies the MSFO option,
      the filename is required on the record which appears first in the
      control  stream.   If filename is present on both, it must be the
      same.  If the filename is  omitted  on  the  second  of  the  two
      records,  the MSF option must be the final option on that record.
      Default: If the MSFO option is not used, no MSF is output.

      If a MSFO is output, then the iteration  estimates  may  also  be
      seen  in the original parameterization for those iterations whose
      summaries appear in intermediate printout.  These  estimates  may
      be found in file INTER.

      When  MAXEVAL=0  and the Covariance Step is implemented, the MSFO
      option may also be used, and then a model specification file will
      be output which will include information from the Covariance Step
      and from the input model specification file concerning  the  ear-
      lier  Estimation  Step (in this case there must be an input model
      specification file).
      (See model_specification_file).

 MUM=s
      s is a string of [M|N|D|X]'s  with  each  letter  representing  a
      THETA  parameter  in  numerical order. M indicates that the THETA
      should be Mu  modeled.  N indicates the THETA should  not  be  Mu
      modeled.  D  indicates  the  program will decide if the parameter
      should be Mu modeled or not.  X indicates that THETA is  involved
      in a covariate-dependent mixture model and is required if this is
      the case.
      Default is DDDD...

 MUM=v1(n1):v2(n2):v3(n3)...
      An alternative syntax may be used.  V is a letter (N,M,D, or  X),
      and  n  is  a number list.  For example, to specify that thetas 3
      and 5 through 8 should not be MU modeled, theta 2 is a population
      mixture parameter, and thetas 6 and 12 are to be MU modeled:
      MUM=N(3,5-8):X(2):M(6,12)

 NBURN=n
      When  used with the SAEM method n is the maximum number of itera-
      tions used to perform the stochastic  phase.   Default  is  2000.
      When used with the BAYES method n is the maximum number of itera-
      tions used to perform the burn-in phase. Default is 4000.

 NITER=n
      When used with the ITS, IMP, and IMPMAP methods n is the  maximum
      number  of  iterations.  Default  is 50.  When used with the SAEM
      method n5 is the number  of  iterations  for  the  non-stochastic
      accumulation  phase.  Default  is  1000. When used with the BAYES
      method n5 is the number of iterations used to obtain the station-
      ary  distribution.   Default  is  10000.  NITER may also be coded
      NSAMPLE.

 NOCOV=[0|1] (NM73)
      If covariance estimation is not desired for a particular  estima-
      tion  step,  set NOCOV=1.  It may be turned on again for the next
      estimation  step  with  NOCOV=0.    If  NOCOV=1  is  set  for  an
      FOCE/Laplace/FO step, this is equivalent to $COV NOFCOV setting.
      For  ITS  and  IMP,  covariance estimation can take some time for
      large problems, and you may wish to  obtain  only  the  objective
      function, such as in the case of $EST METHOD=IMP EONLY=1 after an
      SAEM estimation. NOCOV has no effect on  BAYES  analysis,  as  no
      extra time is required in assessing covariance for BAYES.
      By  default, standard error information for the classical methods
      (FO/FOCE/Laplace) will be given only if they are the last estima-
      tion method, even if NOCOV=0 for an intermediate estimation step.
      If NOCOV=1 for the FOCE/LAPLACE/FO method, and  it  is  the  last
      estimation  step,  then  standard error assessment for it will be
      turned off.

 NOLABEL=[0|1]
      1 indicates that the row of item names for FILE will not be writ-
      ten,  otherwise  0, the default.  Affects the raw output file and
      all additional output files.

 NONINFETA=[0|1] (NM73)
      Sometimes, gradients are not  properly  evaluated  for  classical
      NONMEM methods, when not all etas are used for all subjects.  For
      example, an eta to a ka absorption rate constant during a fit  of
      a  subject with only IV dosing would be considered a non-influen-
      tial data.  If $EST METHOD=1  or  0  is  used  without  the  SLOW
      option,  this  can  result in evaluating very large and incorrect
      gradients, which in turn affects the search path,  and  sometimes
      the  final objective function value.  Should this occur, add NON-
      INFETA=1 to the $EST record.  NONINFETA=0 is the default.

 NOPRIOR=[0|1]
      If prior information was specified using  the  $PRIOR  statement,
      then normally the analysis is set up for three stage hierarchical
      analysis.  By default NOPRIOR=0, and this prior information  will
      be  used.   If NOPRIOR=1, then for the particular estimation, the
      prior information is not included in the analysis. This is useful
      if you do not want to use prior information during a maximization
      (METHOD=IMP, CONDITIONAL, IMPMAP, SAEM, or ITS), but then use  it
      for  the  Bayesian analysis (METHOD=BAYES).  With NOPRIOR=1, FOCE
      is still allowed to evaluate an S MATRIX, since prior information
      is not used.  I.e.,  $EST NOPRIOR=1 and $COV MATRIX=S are permit-
      ted.  With NONMEM 7.3, when NOPRIOR=1 is set, the estimation will
      not  use  TNPRI prior information (TNPRI should only be used with
      FO/FOCE/Laplace estimations). In  previous  versions  of  NONMEM,
      NOPRIOR=1 did not act on TNPRI priors.

 NOSUB=[0|1] (NM74)
      With  NOSUB=0,  label  substitution  will  be performed for final
      estimates in the NONMEM report file.   (See $ABBREVIATED).   This
      is  the  default.   With  NOSUB=1, label substitution will not be
      performed.

 NOTITLE=[0|1]
      1 indicates the header for FILE will not be written, otherwise 0,
      the default.  Affects the raw output file and all additional out-
      put files.

 NUMDER=[0|1|2|3] (NM73)
      With NUMDER=1, NONMEM computes and displays numerically evaluated
      derivatives  of  Y  or  F  with respect to eta and eps (G and H).
      These numeric values are displayed in root.fgh, but are not  used
      in estimation.
      With  NUMDER=0,  file  root.fgh  is  not  produced.   This is the
      default.
      With  NUMDER=2,  analytical  derivatives  values  are  stored  in
      root.agh
      With NUMDER=3, both root.agh and root.fgh are produced.

 NUMERICAL
      Requests  that second eta-derivatives for the Laplacian method be
      obtained numerically.

 NONUMERICAL
      Requests that second eta-derivatives for the Laplacian method  be
      computed  by  PRED.  Not permitted with the combination LAPLACIAN
      and INTERACTION.  Otherwise, this is the default.

 NUTS_BASE (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0.025.

 NUTS_DELTA  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0.8

 NUTS_EPARAM  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0

 NUTS_GAMMA  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0.05

 NUTS_INIT  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0.075.

 NUTS_MASS  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is B

 NUTS_MAXDEPTH  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 10

 NUTS_OPARAM  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 1

 NUTS_REG  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0.0

 NUTS_SPARAM  (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 1

 NUTS_STEPINTER (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0

 NUTS_STEPITER (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 1

 NUTS_TERM (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0.05

 NUTS_TEST (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0

 NUTS_TRANSFORM (NM74)
      See Guide INTRODUCTION TO NONMEM 7. Default is 0

 OACCEPT=n
      Used only with the BAYES method. n has  meaning  only  for  OMEGA
      sampled  by  the  Metropolis-Hastings  algorithm.  The scaling of
      degrees of freedom is adjusted so that  samples  are  accepted  n
      fraction of the time. See OSAMPLE_M2=. Default is 0.5.

 OLKJDF=n (NM74)
      Used  with  NUTS  method,  OLKJDF  stands  for  Omega LKJ density
      degrees of freedom.  When 0, the usual inverse Wishart  prior  is
      used  for  Omegas. When OLKJDF>0, then the LKJ density is used as
      the prior, with OLKJDF degrees of freedom for all  omega  blocks.
      In  addition, only diagonal elements of the OMEGA prior are used,
      assuming a density dependent on the OVARF value.  OLKJDF  may  be
      set  >0  when using METHOD=BAYES as well, but thetas will then be
      M-H sampled using the OSAMPLE_M1, OSAMPLE_M2, and OSAMPLE_M3 set-
      tings.
      Default is 0.
      See  record  $OLKJDF (NM75) to specify LKJ correlation degrees of
      freedom for each omega block.

 OLNTWOPI (NM74)
      The   objective    function    is    reported    including    the
      NETA*NIND*LOG(2pi)  constant  term for SAEM and BAYES, where NETA
      is the number of etas, and NIND is number of  individuals.   Com-
      pare LNTWOPI. Either or both may be used.

 THETABOUNDTEST, OMEGABOUNDTEST, SIGMABOUNDTEST
      With NONMEM VI, the estimation step sometimes terminates with the
      message
      PARAMETER ESTIMATE IS NEAR ITS DEFAULT BOUNDARY.
      These options request that the "default boundary  test"  be  per-
      formed for THETA, OMEGA, and SIGMA, respectively.  THETABOUNDTEST
      may also be coded TBT or TBOUNDTEST; OMEGABOUNDTEST may  also  be
      coded  OBT or OBOUNDTEST; SIGMABOUNDTEST may also be coded SBT or
      SBOUNDTEST.  These options are the defaults.

 NOTHETABOUNDTEST, NOOMEGABOUNDTEST, NOSIGMABOUNDTEST
      Instructs NONMEM to omit the "default  boundary  test"  for  this
      type  of  variable, i.e., to behave like NONMEM V in this regard.
      Any option listed above may be  preceded  by  "NO".   The  THETA,
      OMEGA, and SIGMA choices are independent of each other.  E.g., it
      is possible to specify  NOOBT  (to  prevent  the  "default  OMEGA
      boundary test") and permit both the "default THETA boundary test"
      and "default SIGMA boundary test".

 OMITTED
      The Estimation Step is not implemented.

 OPTMAP=n (NM73)
      For alternative MAP (eta optimization) methods.   With  OPTMAP>0,
      SLOW option may be needed.

      OPTMAP=0
           Standard  variable  metric (Broyden, Fletcher, Goldfarb, and
           Shanno (BFGS)) optimization method used by  NONMEM  to  find
           optimal eta values (referred to as eta hat) for each subject
           at the mode of their posterior densities,  using  analytical
           derivatives  of  F  with respect to etas (G), and analytical
           derivatives of F with respect to etas (H),  that  were  sup-
           plied by NMTRAN or by the user.  This is the default.

      OPTMAP=1
           Variable  metric  method  using  numerical finite difference
           methods for first derivatives of F  with  respect  to  etas.
           Necessary  when  not  all code used in evaluating F, G and H
           for observation event records is abbreviated code (some  may
           be  in  verbatim code), and/or some portions of the computa-
           tion of  F, G and H are evaluated  in  a  hidden  subroutine
           specified by "$SUBROUTINES OTHER=" and the user-written code
           does not compute  the  eta  derivatives.  When  OPTMAP=1  is
           present,  values of G and H are ignored during eta optimiza-
           tion.  This may  be  used  to  test  user-coded  deriatives,
           because  two  runs,  one  with  OPTMAP=1 and one without it,
           should give very similar values for the OBJV, WRES, etc.  if
           the user-coded derivatives are correct.

      OPTMAP=2
           Nelder  Mead method, which uses a secant method, rather than
           relying on derivatives.

 ORDER=xxxf
      The values of x may be T (Theta), S (Sigma), and O (Order).   The
      value of f may be U (Upper) or L (Lower).  Affects the way theta,
      omega, and sigma are displayed in the raw and  additional  output
      files.  xxx gives the overall order, and f gives the order within
      OMEGA and SIGMA.  Affects the raw output file and all  additional
      output  files.  The default is TSOL: THETA, SIGMA, OMEGA in Lower
      triangular form.  Does not affect the NONMEM report file.
      (See order_option).

 OSAMPLE_M1=n
      Used only with the BAYES method. n is the number  of  times  each
      iteration  that  OMEGA is generated using the Metropolis-Hastings
      algorithm and a Wishart proposal density that has variance  based
      on the previous samples. If n < 0 this indicates that the OMEGA's
      are Gibbs sampled using the appropriate Wishart proposal  density
      and other options are not relevant.  Default is -1.

 OSAMPLE_M2=n
      Used  only  with  the  BAYES method. n has meaning only for OMEGA
      sampled by the Metropolis-Hastings algorithm.  n is the number of
      times that OMEGA is generated using a Wishart proposal density at
      the present OMEGA position and degrees of freedom scaled to  have
      samples accepted a particular fraction of the time. If n < 0 this
      is done as many times as  there  are  non-fixed  OMEGA  elements.
      Default is -1.

 OVARF=x (NM74)
      Used  with  NUTS  method  and OLKJDF option.  OVARF is the weight
      factor to STD prior to the log sqrt OMEGA diagonal elements,  the
      normal density of the log square root of OMEGA centered about log
      square root of Omega prior, and scaled with  OVARF  (see  below).
      That            is,            log(sqrt(Omega(i)))           Nor-
      mal(log(sqrt(OmegaPrior(i))),1/OVARF).  If OVARF<0, then a  half-
      t-distribution  of  degrees of ABS(OVARF) is used as the prior to
      the sqrt of OMEGA diagonal elements.  Use OVARF=-1 for the  half-
      Cauchy distribution.
      Default is 1.  See also $OVARF control record.

 PACCEPT=n
      Used  only  with the BAYES method. n has meaning only for popula-
      tion parameters sampled  by  the  Metropolis-Hastings  algorithm.
      The  scaling of variance is adjusted so that samples are accepted
      n fraction of the time. See PSAMPLE_M2=. Default is 0.5.

 PARAFILE=filename
      Name of the "parallel file" (the  parallelization  profile)  that
      controls  parallelization  (distributed computing).  Default file
      name if not specified: parallel.pnm or parafile name specified on
      nmfe command.
      PARAFILE=ON turns on parallelization for this $ESTIMATION record.
      PARAFILE=OFF  turns  off  parallelization  for  this  $ESTIMATION
      record.

 PARAFPRINT=n (NM74)
      The print iteration intervals to the parallelization log file can
      be controlled by this option during parallelization of  the  $EST
      step.   See  also $COVARIANCE record and nmfe74 command.  Default
      is PARAFPRINT=1.

 PHYTYPE=n (NM74)
      Default is 0.  By default, after an estimation is performed,  the
      phi(),  conditional means of the individual parameters, and their
      variances, are reported in the root.phi file, where root  is  the
      root  name of the control stream file. If you wish to have condi-
      tional mean etas reported, set PHITYPE=1.
      See "Stochastic  Approximation  Expectation  Maximization  (SAEM)
      Method" in Guide INTRODUCTION TO NONMEM 7.

 POSTHOC
      This  option  may  be used when the FO method is used.  After the
      Estimation Step terminates, the eta values are estimated for each

      individual.   To estimate the etas based on the initial estimates
      of THETA, OMEGA, and SIGMA (found either in the control stream or
      in  a  model  specification  file), also specify MAXEVAL=0 (which
      omits the Estimation Step).  The  conditional  estimates  of  the
      etas are referred to as Conditional Parametric Etas (CPE).

 NOPOSTHOC
      Etas  are not estimated.  This is the default with METHOD=0.  May
      not be used with METHOD=1.

 PREDICTION
      Indicates that Y (with NM-TRAN abbreviated code)  or  F  (with  a
      user-supplied  PRED  or  ERROR  code)  will serve as a prediction
      variable, i.e., it will be set to a prediction.  Upon simulation,
      the  simulated  observation is possibly also being set in Y or F.
      (However, the DV data item may instead be  set  directly  to  the
      simulated  observation.)  Also, etas (if any) are population etas
      only if epsilons also appear.  This is the default.  Compare with
      LIKELIHOOD, -2LOGLIKELIHOOD options.

 PRINT=n
      Iteration summaries are printed for the 0th, every nth iteration,
      and  last  iteration.   When  n=0,  no  summaries  are   printed.
      Default:  9999  (so  that  summaries are printed for 0th and last
      iterations).

 PRIORC (NM74)
      The objective function is reported including the  prior  constant
      term (constant term to the prior).

 PSAMPLE_M1=n
      Used  only  with the BAYES method. n has meaning only for popula-
      tion parameters sampled by the Metropolis-Hastings algorithm.   n
      is  the  number of times that a vector of THETA's and SIGMA's are
      generated using a multivariate normal proposal density  that  has
      mean/variances based on the previous samples. Default is 1.

 PSAMPLE_M2=n
      Used  only  with the BAYES method. n has meaning only for popula-
      tion parameters sampled by the Metropolis-Hastings algorithm.   n
      is  the  number of times that a vector of THETA's and SIGMA's are
      generated using a multivariate normal proposal density  that  has
      mean  at  the  present  parameter position and variance scaled to
      have samples accepted a particular fraction of the time. If n < 0
      this  is  done  as  many  times  as there are Metropolis-Hastings
      parameters. Default is -1.

 PSAMPLE_M3=n
      Used only with the BAYES method. n has meaning only  for  popula-
      tion  parameters sampled by the Metropolis-Hastings algorithm.  n
      is the number of times in mode 3 that each parameter is individu-
      ally sampled. Default is 1.

 PSCALE_MIN=n, PSCALE_MAX=n (NM74)
      In  MCMC  sampling, the scale factor used to vary the size of the
      variance  of   the   proposal   density   population   parameters
      (theta/sigma)  that  are  not Gibbs sampled, in order to meet the
      PACCEPT condition, is by default bounded by PSCALE_MIN  of  0.01,
      and  PSCALE_MAX=1000.   This should left alone for MCMC sampling,
      but on occasion there may be a reason to  expand  the  boundaries
      (perhaps to PSCALE_MIN=1.0E-06, PSCALE_MAX=1.0E+06).

 RANMETHOD=[n|S|m|P] (NM72)
      See INTRODUCTION TO NONMEM 7, Reference [5] and [7].
      n:  the  random  number generator used for all Monte Carlo EM and
      Bayesian methods.
        0: ran0 of reference [5], minimal standard generator
        1: ran1 of reference [5], Bays and Durham.
        2: ran2 of reference [5].
        3: ran3 of reference [5], Knuth. (Default)
        4: NONMEM's traditional random number generator used in $SIMULATION
      S: sobol method without scrambling,  used  during  importance  or
      direct  sampling  (methods  IMP, IMPMAP, and DIRECT) and only for
      the purpose of creating quasi-random samples of eta vectors.   As
      of  NONMEM  7.3,  Sobol may be used for BAYES and SAEM methods as
      well.
      m: the type of scrambing desired
        0: no scrambing (S0 is the same as S)
        1: Owen type scrambling
        2: Faure-Tezuka type scrambling
        3: Owen plus Faure-Tezuka type scrambling.
      P: each subject will receive its own seed path,  that  will  stay
      with  that subject regardless of whether the job is run as a sin-
      gle process or parallel process. (NM74)
      The RANMETHOD specification propagates to subsequent $EST records
      in  a  given  problem, but does not propagate to $CHAIN or $TABLE
      records.

 REPEAT
      The search is repeated with the initial estimates being the final
      estimates  from  the first search and with new UCP, so that a UCP
      value of 0.1 now corresponds to a final estimate from  the  first
      search.  Cannot be used with STIELTJES.

 NOREPEAT
      The  estimate  obtained  at the end of the minimization search is
      taken to be the final parameter estimate. This  is  the  default.
      Cannot be used with STIELTJES.

 REPEAT1
      The search of the first stage of the Stieltjes method is repeated
      with the initial estimates being the  final  estimates  from  the
      first  search  and  with  new UCP, so that a UCP value of 0.1 now
      corresponds to a final estimate from the first search.  May  only
      be used with STIELTJES.

 NOREPEAT1
      The estimate obtained at the end of the search of the first stage
      of the Stieltjes method is taken to be the final parameter  esti-
      mate  at  the  first stage. This is the default. May only be used
      with STIELTJES.

 REPEAT2
      The search of  the  second  stage  of  the  Stieltjes  method  is
      repeated  with  the  initial  estimates being the final estimates
      from the first search and with new UCP, so that a  UCP  value  of
      0.1  now  corresponds  to a final estimate from the first search.
      May only be used with STIELTJES.

 NOREPEAT2
      The estimate obtained at the end of  the  search  of  the  second
      stage  of the Stieltjes method is taken to be the final parameter
      estimate at the second stage. This is the default.  May  only  be
      used with STIELTJES.

 SADDLE_HESS=n (NM74)
      SADDLE_HESS=0  selects  the  Hessian matrix last generated by the
      variable metric method.  SADDLE_HESS=1 causes the full second de-
      rivative  information  matrix  (identical to R matrix in the $COV
      step) to be evaluated.  Default is 0.
      See "Resetting the Search to  Circumnavigate  Saddle  Points  and
      Detect  Inestimable Parameters" in Guide INTRODUCTION TO NONMEM 7
      for a discussion of SADDLE_HESS and SADDLE_RESET options.

 SADDLE_RESET=n (NM74)
      SADDLE_RESET is the number of times that a reset should occur  in
      the course of the search.  Normally, should be set to 1.  Default
      is 0.

 SEED=n
      The initial seed for the random number  generator  used  for  the
      Monte-Carlo methods. Default is 14455.

 CLOCKSEED=[0|1] (NM75)
      If  CLOCKSEED=1  (default  is  0),  actual  starting seed will be
      10000*(seconds  after  midnight)+SEED.   This  allows  a  control
      stream  to  produce  different  stochastic  results for automated
      replications, without the need to modify the seed  value  in  the
      control stream file in each replication.

 SELECT=[0|1|2|3] (NM73)
      Used  with  METHOD=CHAIN  and $CHAIN to specify how the sample is
      selected.
      SELECT=0
      If ISAMPEND>=ISAMPLE,  then  the  default  action  for  selecting
      between   ISAMPLE   and   ISAMPEND   is  taken,  which  for  $EST
      METHOD=CHAIN is to find the one giving the best OBJ at  the  ini-
      tial  values, and for $CHAIN is to randomly select a sample, with
      replacement.  This is the default.
      SELECT=1
      The sample is selected sequentially from ISAMPLE to ISAMPEND with
      each  new  use of $CHAIN/$SIML with multiple sub-problems for the
      given problem, and with each new $EST METHOD=CHAIN with  multiple
      sub-problems  and  across problems. When ISAMPEND is reached, the
      sample selection begins at ISAMPLE again.
      SELECT=2
      Uniform random selection of sample, without replacement.   Should
      the sample selection become exhausted, which would occur if CHAIN
      or $CHAIN records are utilized for more  than  ISAMPEND-ISAMPLE+1
      times, subsequent sample selection then occurs with replacement.
      SELECT=3
      Uniform  random  selection  of  sample, with replacement (this is
      equivalent to SELECT=0 for $CHAIN).

 SIGDIGITS=n
      Number of significant digits  required  in  the  final  parameter
      estimate.   SIGDIGITS  is  not  used  by the Monte-Carlo methods.
      Default: 3.  May also be coded NSIGDIGITS.

 SLKJDF=n (NM74)
      Used with NUTS  method,  SLKJDF  stands  for  Sigma  LKJ  density
      degrees  of  freedom.  When 0, the usual inverse Wishart prior is
      used for Sigmas. When SLKJDF>0, then the LKJ density is  used  as
      the  prior,  with  SLKJDF  degrees of freedom.  In addition, only
      diagonal elements of the Sigma prior are used.  SLKJDF may be set
      >0  when using METHOD=BAYES as well, but Sigmas (in cholesky for-
      mat) will then be M-H sampled using the  PSAMPLE_M1,  PSAMPLE_M2,
      and  PSAMPLE_M3 settings (choleskys of sigma elements are treated
      as extensions of the THETA parameters in M-H sampling methods).
      Default is 0.
      See record $SLKJDF (NM75) to specify LKJ correlation  degrees  of
      freedom for each sigma block.

 SIGL=n
      n is used to calculate the step-size for finite difference deriv-
      atives independent of the SIGDIGITS value.  If n=0 or n=100  then
      SIGL  is  ignored  and  SIGDIGITS is used as in versions prior to
      NONMEM 7.  SIGL should usually be 2 to 3 times the value of NSIG.
      It is not used by the SAEM or BAYES methods.

 SIGLO=n (NM72)
      The  precision  to  which the individual etas are optimized.  The
      SIGL value set by the user continues  to  be  the  precision  (or
      delta  )  setting  for  the  finite  difference algorithms in the
      higher level estimation process for THETAS, OMEGAS,  and  SIGMAS.
      By  default,  if SIGLO is not specified, then SIGLO is set to the
      same value as SIGL.  Should SIGLO be used, the  recommended  set-
      ting would be:
      SIGLO<=TOL
      SIGL<=SIGLO
      NSIG<=SIGL/3

 SIGMABOUNDTEST
      See OMEGABOUNDTEST.

 SLOW Requests  a slower method of computation.  Required when either a
      mixture model is used along with CENTERING, or NUMERICAL is used.
      If not present, the option is automatically supplied in these two
      cases.  For problems where NONMEM VI  does  not  behave  as  well
      (e.g.  yields  a higher OFV at termination) compared to NONMEM V,
      inclusion of the SLOW option may sometimes, but not always, yield
      NONMEM VI results that are similar to NONMEM V.

 SLOW=1
      Same as SLOW.

 NOSLOW
      Requests  a  faster  method  of computation.  This is the default
      (but see SLOW)

 SLOW=2
      This option is permitted with STIELTJES.

 SORT Individual contribution to the objective function value and indi-
      vidual  contributions to the gradients are sorted before they are
      summed, so that smaller numbers are summed before larger numbers.

 NOSORT
      Individual contribution to the objective function value and indi-
      vidual  contributions to the gradients are summed in the order in
      which the individual records appear in the NONMEM  data  set,  as
      was done prior to NONMEM VI.  This is the default.

 STDOBJ=x (NM74)
      For  importance sampling and direct sampling only, if ISAMPEND is
      specified as an upper integer value, and STDOBJ is set to a  real
      value  greater  than 0, then NONMEM will vary the number of Monte
      Carlo samples under each subject between  ISAMPLE  and  ISAMPEND,
      until the stochastic standard deviation of the objective function
      falls below STDOBJ.

 STIELTJES
      A set of tentative population estimates are first obtained  using
      some  1st-  or 2nd-order method.  A tentative value for the inte-
      gral (i.e. an area) is obtained.  Then numerical  integration  is
      used  to  obtain second-stage estimates.  See the Introduction to
      Version VI 2.0.  Not permitted with METHOD=HYBRID.

 SVARF=x (NM74)
      Used with NUTS method and SLKJDF option.   SVARF  is  the  weight
      factor  to STD prior to the log sqrt SIGMA diagonal elements, the
      normal density of the log square root of SIGMA centered about log
      square  root  of  SIGMA prior, and scaled with SVARF (see below).
      That           is,            log(sqrt(Sigma(i)))            Nor-
      mal(log(sqrt(SigmaPrior(i))),1/SVARF).   If SVARF<0, then a half-
      t-distribution of degrees of ABS(SVARF) is used as the  prior  to
      the  sqrt of SIGMA diagonal elements.  Use SVARF=-1 for the half-
      Cauchy distribution.
      Default is 1.
      See also $SVARF control record.

 TBLN=n (NM75)
      Used with $EST METHOD=CHAIN and $CHAIN records to allow selecting
      a  table  within a raw output file. See "Method for creating sev-
      eral instances for a problem  starting  at  different  randomized
      initial positions" in Guide INTRODUCTION TO NONMEM 7.

 THETABOUNDTEST
      See OMEGABOUNDTEST.

 THIN=n (NM74)
      The  Bayesian  records  retained  in  the  raw output file may be
      adjusted by every THINth iteration. So, if  THIN=10,  then  every
      10th  iteration  is  recorded  in  the raw output file. The PRINT
      option controls only the iterations printed to  the  console  and
      NONMEM report file.
      Default is 1.

 TPU=n (NM75)
      If TPU>0, use THETA_PRIORU routine in ..\source\THETA_PRIORU.f90.
      You can make a copy of THETA_PRIORU.f90, modify it, call it USER-
      PRIORT.f90,  for  example,  then  specify it as an OTHER routine:
      $SUBR OTHER=USERPRIORT.f90.  See also \$OLKJDF and  \$SLKJDF  for
      OMEGA and SIMGA.
      Default is 0.

 TTDF=n (NM74)
      TTDF  stands  for Theta t-density degrees of freedom.  It is used
      with NUTS method.  Default is 0.  When 0, the usual  normal  den-
      sity prior is used as a prior density for thetas.  When TTDF>0, a
      t-distributed prior is used.  TTDF  may  be  set  >0  when  using
      METHOD=BAYES  as  well, but thetas will then be M-H sampled using
      the PSAMPLE_M1, PSAMPLE_M2, and PSAMPLE_M3 settings.  TTDF may be
      a real number.
      See  also $TTDF control record (NM75) to specify degrees of free-
      dom for each theta.  The value of TTDF overrides $TTDF.

 ZERO=list
      Required with METHOD=HYBRID.  A list of indices  for  etas  which
      are  fixed  to  zero during the Estimation Step.  "list" contains
      one or more integers.  If more than one, they must be  surrounded
      by parentheses.  The list must be contained entirely on one line.
      The indices may be separated by commas or spaces.

 Reserved Variables that are of Interest During the Estimation Step

 MUFIRSTREC (NM74)
      The MUFIRSTREC option can speed up  the  NUTS  method,  and  also
      ordinary  BAYES,  FAST  FOCE,  ITS,  and  the  EM  methods.   Set
      MUFIRSTREC=1 in $PRED or $PK. MUFIRSTREC=1 selects the  covariate
      of  the  first record of the subject, rather than averaging among
      its records when using that covariate in a MU reference.
      The first statement of the $PRED or $PK block should be
      include nonmem_reserved_general.

 OBJQUICK (NM74)
      The OBJQUICK option can speed up the NUTS method, and also  ordi-
      nary  BAYES, FAST FOCE, ITS, and the EM methods.  Set OBJQUICK in
      $PRED or $PK.

      OBJQUICK=0
           Default. Standard NONMEN processing of the model occurs.

      OBJQUICK=1
           Certain tests and initializations are skipped.

      OBJQUICK=2
           A simplified modeling process occurs, but  which  cannot  be
           used when $LEVEL or $MIX is used in the model.
      The first statement of the $PRED or $PK block should be
      include nonmem_reserved_general.

 REFERENCES: Guide I, section C.3.5.1 
 REFERENCES: Guide II, section E , F 
 REFERENCES: Guide IV, section III.B.14 , IV.G 
 REFERENCES: Guide V, section 9.4.1, , 13.2 
 REFERENCES: Guide VII, section I , II , III 
 REFERENCES: Guide Introduction_7
