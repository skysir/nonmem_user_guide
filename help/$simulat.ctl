


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            $SIMULATION                             |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Instructions for the NONMEM Simulation Step
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SIMULATION  (seed1 [seed2] [NORMAL|UNIFORM|NONPARAMETRIC] [NEW]) ...
              [CLOCKSEED=[0|1]]
              [SOURCE_EPS=n]
              [SUBPROBLEMS=n] [ONLYSIMULATION] [OMITTED]
              [REQUESTFIRST] [REQUESTSECOND] [PREDICTION|NOPREDICTION]
              [TRUE=INITIAL|FINAL|PRIOR] [TTDF=n]
              [BOOTSTRAP=n [REPLACE|NOREPLACE] [STRAT=label] [STRATF=label]]
              [NOREWIND|REWIND] [SUPRESET|NOSUPRESET]
              [RANMETHOD=[n|S|m|P] ]
              [PARAFILE=[filename|ON|OFF]

 SAMPLE:
 $SIMULATION     (889215690) (2239177789 UNIFORM)

 DISCUSSION:
 Optional.  Requests  that  the  NONMEM Simulation Step be implemented.
 May also be coded $SIMULATE or $SIML.

 Usually, when the Simulation Step is implemented, the simulated obser-
 vation  is  taken to be the quantity to which the Y variable (with NM-
 TRAN abbreviated code) or F variable (with  a  user-supplied  PRED  or
 ERROR  routine)  is set.  This is the default behaviour.  However, see
 option NOPREDICTION below.

 If a $ESTIM record appears in the problem specification,  then  unless
 the LIKELIHOOD or -2LOGLIKELIHOOD option appears on the $ESTIM record,
 etas (if any) are understood to be single-subject  etas,  except  when
 epsilons also appear, in which case the etas are understood to be pop-
 ulation etas.  If a $ESTIM record does not appear, but a $SIMUL record
 appears,  then  unless  the  NOPREDICTION option appears on the $SIMUL
 record, etas are understood in the same way.   When  the  NOPREDICTION
 option is used, the etas are understood to be population etas.

 In NM-TRAN abbreviated code, there can be a special block of code that
 is implemented only during the simulation task  (See ICALL, Simulation
 Block).  It is called a "simulation block".

 OPTIONS:
 The  information  coded  within  each  set  of parentheses defines the
 attributes of a single random source.  A source of random  numbers  is
 an  "infinite"  stream of random numbers.  Each pair of parentheses on
 the record defines a separate source of random numbers, and the infor-
 mation  coded  within  the  parentheses  defines the attributes of the
 source.  The sources are ordered as they are  defined  on  the  $SIMUL
 record.   The  numbers  from  a source are explicitly available to the
 user via the NONMEM utility routine: RANDOM (See RANDOM).  By default,
 the first source is earmarked for the simulation of etas and epsilons,
 and then the numbers from this source are not explicitly available  to
 the user.  However, see option SOURCE_EPS below.

 seed1
      Seed1 is the first seed for the random source, an integer between
      0 and 2147483647.  If this is not the first problem specification
      in  the control stream, then seed1 can be -1, indicating that the
      random source is to be continued from the previous problem.

 seed2
      Seed2 is the second  seed  for  the  random  source,  an  integer
      between  0  and 2147483647.  For use of a second seed, see NONMEM
      Users Guide, Part IV.

 NORMAL
      The random numbers of the source are  to  be  pseudo-normal  with
      mean 0 and variance 1 (unless the source is the first and used to
      generate eta and epsilon realizations, in which  case  the  vari-
      ance-covariance  of  these  variables  is  that  specified in the
      $OMEGA and $SIGMA records).  This is the default.

 UNIFORM
      The random numbers of the source are to be pseudo-uniform on  the
      interval [0,1].

 CLOCKSEED=[0|1]  (NM75)
      If  CLOCKSEED=1  (default  is  0),  actual  starting seed will be
      10000*(seconds after midnight)+SEED (applies to  both  seed1  and
      seed2,  if  specified.)   This allows a control stream to produce
      different stochastic results for automated replications,  without
      the  need  to modify the seed value in the control stream file in
      each replication.

 SOURCE_EPS=n (NM75)
      May be used to specify that the simulation of epsilons should use
      a different source than the default, which is the first. n is the
      number of the source for epsilons.  This  source  must  have  the
      NORMAL attribute.

 NONPARAMETRIC
      The  random  numbers  from the first source defined with the NON-
      PARAMETRIC attribute are used to generate realizations of  random
      vectors  from a (multivariate) nonparametric distribution on eta,
      obtained from the Nonparametric Step of an earlier  problem.   An
      input  MSF  must also be present.  It is understood that the etas
      are to be simulated from the  nonparametric  distribution  rather
      than  from  the  pseudo-normal  distribution  associated with the
      first source.  The NONPARAMETRIC attribute can only  be  used  in
      the definition of the second or subsequent source.

 NEW  If  the  NEW  option  is  used,  the  vector of eta's (epsilon's)
      changes with each call to SIMETA (SIMEPS), rather  than  only  at
      the  start of the next individual record (next level-two record).
      That is, with NEW, a call to SIMETA will obtain new  eta's,  even
      if  ID  has not changed.  A call to SIMEPS will obtain new eps's,
      even if L2 has not changed.  The output NONMEM  report  describes
      this  as  "DIFFERENT  ETA  AND  EPS  WITH EACH CALL TO SIMETA AND
      SIMEPS".

 ONLYSIMULATION
      NONMEM is being asked to simulate data but  not  to  evaluate  an
      objective  function  on  these  data.   WRES values in tables and
      scatterplots will be 0.  PRED-defined data items  in  tables  and
      scatterplots  will  be  computed using simulated etas and initial
      thetas.

      $ESTIM, $COV and $NONP cannot be used with ONLYSIMULATION.  Also,
      see the PREDICTION and NOPREDICTION options.

 SUBPROBLEMS=n
      Requests that the entire NONMEM problem is to be repeated n times
      in succession (including all NONMEM  steps:  simulation,  estima-
      tion,  covariance, table, scatterplot).  Each subproblem includes
      the Simulation Step, but the random sources are simply  continued
      from  subproblem to subproblem.  If n=0 or n=1, there is only one
      subproblem; this is the default.  May  also  be  coded  SUBPROBS,
      NSUBPROBLEMS,  NSUBPROBS.   With all versions of NONMEM, the data
      set for each subproblem after the first is the same data set used
      by the previous subproblem, and includes any changes (transgener-
      ation) made by the previous subproblem.
      With NONMEM 7.4 and higher, see REWIND, below.                    |
      With NONMEM 7.3 and higher, the maximum number of subproblems  is
      increased from 9999 to 2147483647.

 REQUESTFIRST
      NONMEM sets a variable IFIRSTEM in Module ROCM_INT (referenced as
      FIRSTEM in abbreviated code) informing PRED whether or  not  PRED
      needs  to  compute first-partial derivatives with respect to eta.
      Normally, during the Simulation Step, these derivatives  are  not
      needed,  either  by NONMEM or by the user.  However, the user may
      want the first-partial eta derivatives of a PRED-defined item and
      may  want FIRSTEM to reflect this.  With the REQUESTFIRST option,
      the FIRSTEM variable is set so to inform PRED  that  the  deriva-
      tives  need to be computed.  In this case, if an abbreviated code
      is used to compute the PRED-defined item, the item should not  be
      computed within a simulation block, because NM-TRAN does not pro-
      vide derivatives for PRED-defined items in a simulation block.

 REQUESTSECOND
      NONMEM sets a variable ISECDER in Module ROCM_INT (referenced  as
      MSEC  in  abbreviated  code)  informing  PRED whether or not PRED
      needs to compute second-partial derivatives with respect to  eta.
      Normally,  during  the Simulation Step, these derivatives are not
      needed, either by NONMEM or by the user.  However, the  user  may
      want  the  second-partial  eta derivatives of a PRED-defined item
      and may want  the  MSEC  variable  to  reflect  this.   With  the
      REQUESTSECOND  option, the MSEC variable is set so to inform PRED
      that the derivatives need to be computed.  In this  case,  if  an
      abbreviated  code  is  used to compute the PRED-defined item, the
      item should not be computed within a  simulation  block,  because
      NM-TRAN  does not provide derivatives for PRED-defined items in a
      simulation block.  REQUESTSECOND implies REQUESTFIRST.

 PREDICTION
      Permitted only with ONLYSIM, and is the default.
      With or without ONLYSIM, unless the  NOPREDICTION  is  used,  the
      simulated  observation is taken to be the quantity to which the Y
      variable (with NM-TRAN abbreviated code) or F  variable  (with  a
      user-supplied  PRED  or  ERROR  routine) is set.  In a simulation
      block, the DV variable may  be  directly  set  to  the  simulated
      observation, but the Y (or F) variable should also be set to this
      observation.  E.g., if a line of code DV=... is used in a simula-
      tion  block, be sure to follow this line with the additional line
      Y=DV.

 NOPREDICTION
      Permitted only with ONLYSIM.
      Indicates that the simulated observation will be taken to be  the
      value to which the DV variable is set.  The code Y=... is permit-
      ted inside or outside  a  simulation  block,  but  if  such  code
      appears in a simulation block, be sure to also include e.g. DV=Y.
      Also, etas (if any) are understood to be population etas, even if
      epsilons do not appear.

 TRUE=INITIAL
      The initial estimates given in the control stream are used as the
      parameter values ("true values") in the simulation, except when a
      model  specification  file  is input.  When a model specification
      file is input, the initial estimates with the  previous  run  are
      used  as  the parameter values ("true values") in the simulation,
      and the final estimates with the previous run  are  used  as  the
      initial  estimates  in  all tasks other than the simulation.  The
      UCP used with these other tasks are the same as with the previous
      run.   This  is  the default.  May not be used with $MSFI in con-
      junction with SUBPROBLEMS=n (n > 1).

 TRUE=FINAL
      An input model specification file must be used.  The final  esti-
      mates  with  the  previous  run  are used as the parameter values
      ("true values") in the simulation and as the initial estimates in
      all  tasks  other  than  the simulation.  The UCP used with these
      other tasks are new, so that a UCP value of 0.1  now  corresponds
      to a final estimate from the previous run.

 TRUE=PRIOR
      The  values stored in THET_P, OMEG_P, SIGM_P by the PRIOR routine
      are used as the true parameter values ("true values") in the sim-
      ulation.

 TTDF=n (NM75)
      If  TTDF  is  set to a non-zero integer value n, thetas are simu-
      lated with n degrees of freedom t-distribution.  TRUE=PRIOR  must
      be  specified  as  well.  Priors must also be specified appropri-
      ately.
      See also $TTDF control record (NM75) to specify degrees of  free-
      dom for each theta.  The value of TTDF overrides $TTDF
      See INTRODUCTION  TO NONMEM 7, Simulating THETAS with t-Distribu-
      tion

 BOOTSTRAP=n
      With the BOOTSTRAP option, NONMEM does not perform the usual sim-
      ulation  activity  of  randomly creating DV values for a new data
      set, but rather selects a random set of subjects from an existing
      "template"  data  set (which must already have legitimate DV val-
      ues).  The BOOTSTRAP number n refers to how many subjects are  to
      be randomly selected from the data set.  Setting -1 means to ran-
      domly select as many subjects as are in the data set.  For  exam-
      ple,  if  400  subjects  are in the simulation template data set,
      then 400 subjects are randomly selected.  The random  source  is,
      in  effect, uniform, because any subject is equally probable.  If
      n is greater than the number of subjects,  NONMEM  will  use  the
      number of subjects.

      The  BOOTSTRAP option in $SIML is most suitable for multi-subject |
      data, in which there is an ID data column  identifying  the  sub- |
      jects.   However,  see  the example "BOOTSTRAPPING SINGLE SUBJECT |
      DATA" in the INTRODUCTION TO NONMEM 7.

 REPLACE
      Subjects are selected with replacement.   This  results  in  some
      subjects  not  being  selected at all, and some subjects selected
      more than once. This is the default.

 NOREPLACE
      Subjects are  selected  without  replacement,  that  is,  without
      repeating  a  subject.   The  NOREPLACE  feature is reasonable if
      there are many more than n subjects to choose from  in  the  tem-
      plate  dataset  (for  example, 1000 subjects in the template, and
      for each sub-problem, 50 of  them  are  randomly  chosen  without
      replacement, that is, without repeating a subject).

 STRAT
      The label of a data item that serves as the stratification.  This
      splits the data set into distinct sub-sets, guaranteeing  a  spe-
      cific number of subjects will be selected from each category.

 STRATF
      The  label  of a data item that contains the fraction that should
      represent a category  in  the  bootstrapped  data  set.   Without
      STRATF,  the number of subjects to be taken from a given category
      is proportional to the number of subjects in the base data set.

 NOREWIND|REWIND
      The NOREWIND option requests that if any input data  item(s)  are |
      changed  (transgenerated)  during a the simulation step of a sub- |
      problem, they remain changed at the start of the  next  sub-prob- |
      lem.  This  is  the  default.  With NONMEM 7.4, the REWIND option |
      requests that the original data set be used for all sub-problems. |
      If  the  data  set is not changed during simulation, the NOREWIND |
      and REWIND options give the same results.                         |

      Keep in mind that transgeneration performed on the data set using |
      $INFN  with  ICALL=1  affects  the  original  data  set and these |
      changes are unaffected by   REWIND  and  NOREWIND  options.   For |
      example:                                                          |
      $INFN                                                             |
      IF (ICALL==1) THEN                                                |
      DOWHILE(DATA)                                                     |
       ... data transgeneration statements here                         |
      ENDDO                                                             |
      ENDIF                                                             |

 SUPRESET|NOSUPRESET
      This  option  affects seeds for random sources when a $SIMULATION |
      record is included in the scope of a $SUPER problem.  The  SUPRE- |
      SET  option requests that, with subsequent iterations of a super- |
      problem, the seed(s) for all random sources  are  reset  back  to |
      that listed in the $SIMULATION record of the control stream file. |
      This is the default.  With NONMEM 7.4, the NOSUPRESET option  may |
      be used to prevent resetting the random seeds, so that each iter- |
      ation serves as a new random instance.

 RANMETHOD=[n|S|m|P]
      As of NONMEM 7.3, the RANMETHOD option is available for the  $SIM
      record,  to  use alternative random numbers generators.  N is the
      random number generator type, S is Sobol sequence, and m  is  the
      Sobol  scrambler.   The default is NONMEM's traditional one, n=4.
      Among the Sobol sequence methods, the S2 method appears  to  pro-
      vide the least biased random samples, that is nearly uniform dis-
      tribution, with good mixing in multi-dimensional spaces.   As  of |
      NONMEM  7.4,  RANMETHOD  will  also act on the P modifier , which |
      will retain separate seed sequences for each subject, so that the |
      random  variable  patterns are retained regardless of whether the |
      simulation is done in  single  computing  or  parallel  computing |
      mode.
      (See $ESTIMATION).

 PARAFILE=filename
      As  of  NONMEM 7.4, the Simulation Step computation may be paral-
      lelized.  by default. By default, parallelization is  not  turned
      on,  because  simulation is very rapid anyway, and often does not
      need to be accelerated.  When modeling with super-ID  nested  ETA
      levels  ($LEVEL  record  is  present),  parallelization  will not
      occur, since these etas are shared across individuals, and  there
      is  no  guarantee that all subjects sharing the same etas will be
      simulated by the same process.
      PARAFILE=filename specifies a different parafile  than  was  used
      for the previous step.
      PARAFILE=ON turns on parallelization for the Simulation Step.
      Set RANMETHOD=P to permit assure consistient seed patterns
      regardless of whether or not parallelization is performed.  E.g.,
      $SIML ... PARAFILE=ON RANMETHOD=P
      PARAFILE=OFF turns off parallelization for the Simulation
      Step.  This is the default.

 OMITTED
      The Simulation Step is not implemented.

 REFERENCES: Guide IV, section III.B.13 , IV.I 
 REFERENCES: Guide V, section 12.4.8 
 REFERENCES:  Guide  VI, section III.C , III.E , IV.B , IV.G.1 
