


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           $DESIGN (NM75)                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING:  Instructions  for Clinical Trial Design Evaluation and Opti-
 mization
 CONTEXT: NM-TRAN Control Record

 DISCUSSION:
 The optimal design process can help you in designing or  evaluating  a
 clinical  trial.  It may be desired to evaluate specified time points,
 or find the optimal time points,  dose  levels,  and  number  of  time
 points  appropriate for a particular sub-design, etc. The design algo-
 rithms have been modeled after POPED by Hooker et  al.,  and  PFIM  by
 Mentre et al. (see references [22-27]) in Introduction to NONMEM 7.

 USAGE:

 $DESIGN  Record Options
         [APPROX=[FO|FOI|FOCE|LAPLACE|LAPLACEI]]
         [OFVTYPE=[0|1|2|3|4|5|6|7|8]]
         [GROUPSIZE=n]
         [FIMTYPE=[0|1|2|3]]
         [FIMDIAG=[0|1|2|3]]
         [VARCROSS=[0|1]]
         [EOPTD=1]
         [SEED=n]
         [CLOCKSEED=[0|1]]
         [MODE=[0|1|2|3]]
         [DATASIM=[0|1]]
 Additional Control Options for $DESIGN
         [SIGL=n]
         [SIGLO=n]
         [ABORT|NOABORT|NOHABORT]
         [MAXEVALS=n]
         [PRINT=n]
         [NUMERICAL|NONUMERICAL]
         [LIKELIHOOD|-2LOGLIKELIHOOD|-2LL]
         [SLOW|NOSLOW|FAST]
         [POSTHOC|NOPOSTHOC]
         [NOPRIOR=[0|1]]
         [FILE=filename]
         [FORMAT|DELIM=s3]
         [FNLETA=n]
 Options for Setting up Types of Optimal Design
         [NELDER]
         [FEDOROV]
         [RS]
         [STGR]
         [DISCRETE]
         [DISCRETE_RS]
         [DISCRETE_SG]
         [DESEL=label]
         [DESELSTRAT=label]
         [DESELMIN=label]
         [DESELMAX=label]
         [NMIN=label]
         [NMAX=label]
         [STRAT=label]
         [STRATF=label]

 SAMPLE:

 $DESIGN APPROX=FOCEI MODE=1 NELDER FIMDIAG=0
         DATASIM=1 GROUPSIZE=32 OFVTYPE=0

 will  produce  the  most  empirical, "clinical trial simulation" (CTS)
 style covariances, complete with simulated etas and eps, and  standard
 FIM  is  assessed.  If FIMDIAG>0, then a y-expectation covariance will
 be evaluated, but mode will be evaluated with the simulated data.

 DISCUSSION:

 The record name can be shortened  to  $DESI.   Another  name  for  the
 record is $OPTDESIGN, which can be shortened to $OPT.

 The  following is a brief summary of the options.  See Introduction to
 NONMEM 7 "Clinical Trial Design Evaluation and Optimization" for  more
 detail  and more examples.  The options are described in the following
 sections:

 (1) $DESIGN Record Options
 (2) Additional Control Options for $DESIGN
 (3) Options for Setting up Types of Optimal Design

 (1) $DESIGN Record Options

      APPROX=FO (default)
           The nlme approximation method is specified here. First order
           (no interaction) is the default, and the appropriate type of
           covariance matrix is evaluated or used  in  optimization  of
           the design.  Other options are:

           FOI  First order interaction with interaction"

           FOCE first order conditional estimation

           FOCEI
                first order conditional estimation with interaction

           LAPLACE
                Laplace conditional estimation

           LAPLACEI
                Laplace conditional estimation with interaction

      OFVTYPE=1 (default)
           The  objective  function  types  are  comparable to those of
           PopED:

                0,1,3,4,5: design type: d-optimality, -log(det(FIM)), where FIM=Fisher
                Information Matrix (inverse of variance-covariance).
                2: design type: a-optimality, -1/tr(1/FIM)
                6: design type: ds-optimality, -log(det(FIM))+log(det(FIMunintersting))
                To identify parameters as uninteresting, place UNINT at the parameter
                in the same manner you would place FIX.
                7: design type: r-optimality (relative standard error),-1/tr(sqrt(1/FIM)/Paramer))
                8: optimal design type: Individual Bayesian FIM, -log(det(Bayes
                FIM)), as described in the PFIM 4.0 manual [26].

                9: Same as option 2, and using the UNINT filter.
                10: Same as option 7, and using the UNINT filter.

      GROUPSIZE=1 (default)
           The GROUPSIZE is comparable to that of POPED, in  which  the
           FIM is multiplied by this number to provide the subject num-
           ber size of the dataset template.  For  a  template  of  one
           subject,  GROUPSIZE would then offer the variance-covariance
           expected from GROUPSIZE number of subjects.

      FIMTYPE or FIMDIAG=0
           FIMTYPE or FIMDIAG may be set to the following,  and  corre-
           sponds to fimcalc.type in POPED:

           0 (default):  Full information matrix
                Uses finite difference assessment for Theta, Omega, and
                Sigma variances and covariances.

           1  Create a block diagonal information matrix of the esti-
                mates

           2  Create a block information matrix

           3  Create a full information matrix

      VARCROSS=0 (default)
           Standard NONMEM residual variance modeling.

      VARCROSS=1
           Residual Standard Deviation Modeling
           VARCROSS=1 means to treat the residual variance model in the
           manner of PFIM 4.0, as described in the manual.

           FIMTYPE=1 VARCROSS=1 is  equivalent  to  fim.calc.type=4  in
           POPED, and diagonal option in PFIM.

      EOPTD=1
           For  each iteration, this creates a random sample of thetas,
           omegas, or sigmas,  using  the  prior  information.   $PRIOR
           NWPRI  prior information is required, and PLEV=0.999 must be
           specified.  Best used with STGR. See example optdesign16.

      SEED=223345
           Set the starting seed for any random samples to be  created,
           whether for EOPTD=1, or for FOCE type FIM in creating random
           etas (see below).

      CLOCKSEED=0 (default)
           As of nm75, the actual starting seed will be  10000*(seconds
           after  midnight)+SEED  (SEED  may  be set to 0 for this pur-
           pose), if CLOCKSEED=1.  This allows a control stream to pro-
           duce  different  stochastic  results  for automated replica-
           tions, without the need to modify the seed value in the con-
           trol stream file in each replication.

      MODE=0 (default), 1, 2, or 3
           Used  for  specifying  data  and  prediction value type when
           specifying APPROX=FOCEI.  In NONMEM report:
           0 EVALUATE AT F(ETAsim) DURING CONDITIONAL DESIGN ASSESSMENT
           1 EVALUATE AT THE MODE F(ETAhat) DURING CONDITIONAL DESIGN ASSESSMENT
           2 EVALUATE AT F(ETAsim)-G*ETAsim DURING CONDITIONAL DESIGN ASSESSMENT
           3 EVALUATE AT F(0) DURING CONDITIONAL DESIGN ASSESSMENT
           In other words,
           MODE=0 means FOCEI with data  at  F(ETAsim),  and  predicted
           function  evalauted at f(ETAsim), is to be used. This method
           works well.
           MODE=1 means FOCEI with data  at  F(ETAsim),  and  predicted
           function  evaluated  at  the mode, f(ETAhat), is to be used.
           The results are not satisfactory.
           MODE=2   means    linearized    FOCEI,    with    data    at
           F(ETAsim)-G*ETAsim,  and  predicted  function  at F(ETAsim).
           Works well.

      DATASIM=0 (default)
           Normally, y-expectation evaluation of the FIM is  performed.
           To  actually  simulate  data,  set DATASIM=1, and along with
           APPROX=FOCEI, this will produce simulated etas as well.

 (2) Additional Control Options for $DESIGN

      The following options may be set within the $DESIGN  record,  and
      they  operate  exactly  as  their equivalents in the $EST record.
      Thus, the following are equivalent:

      $DESIGN FIMDIAG=1 OFVTYPE=6 APPROX=FO
              MAXEVAL=9999 NOHABORT PRINT=20 SIGL=10 POSTHOC

      $DESIGN FIMDIAG=1 OFVTYPE=6 APPROX=FO
      $EST    MAXEVAL=9999 NOHABORT PRINT=20 SIGL=10 POSTHOC

      The options are as follows.

      SIGL
      SIGLO
      SIGL
      SIGLO
      ABORT/NOABORT/NOHABORT
      MAXEVAL
      MAXEVAL=0 indicates design evaluation (the default)
      MAXEVAL>0 indicates design optimization
      PRINT  control iteration printing during optimal design
      NUMERICAL/NONUMERICAL
      -2LL/LIKELIHOOD/LLIKELIHOOD
      SLOW/NOSLOW/FAST
      POSTHOC
      NOPRIOR
      FORMAT
      FILE
      FNLETA

 (3) Options for Setting up Types of Optimal Design

      The additional options for $DESIGN listed below are for  optimiz-
      ing  parts  of  the  design  components.  For example, the DESEL,
      DESELSTRAT, DESELMIN, DESELMAX can be specified for all the vari-
      ous  design  elements  that you want optimized.  It might be TIME
      for time samples, AMT for dose, or some type  of  covariate  spe-
      cific for the problem.   Certainly, any combination of covariates
      can be requested to be optimized.

      NELDER
           Use Nelder method to search for optimal  continuous  parame-
           ters

      FEDOROV
           Use  to find ideal set of discrete time points from a larger
           set of possible time points.

      RS   Random Search method to find optimal continuous parameters

      STGR Stochastic Gradient method to find optimal continuous param-
           eters

      DISCRETE
           Find optimal number of time points for each sub-design (sub-
           ject template), and use NELDER method to find  optimal  con-
           tinuous parameters.

      DISCRETE_RS
           Find optimal number of time points for each sub-design (sub-
           ject template), and use RS method to find optimal continuous
           parameters.

      DISCRETE_SG
           Find optimal number of time points for each sub-design (sub-
           ject template), and use STGR method to find optimal continu-
           ous parameters.

      Specific parameters must be specified to be optimized.

      This is done using the following options:

      DESEL=label
           The  data  item  (column)  that  contains the design element
           (DESEL) values that are to be modified and  optimized.   For
           example,  TIME  column would indicate that you want times to
           be estimated.

      DESELSTRAT=label
           The  data  item  (column)  indicating  stratification.   The
           DESELSTRAT data item should contain integer indices to indi-
           cate what values are to be shared, and  estimated  together.
           If  a record contains a value of 0 in the DESELSTRAT column,
           then this record is not included in the estimation  process,
           and its value (say its time value in DESEL=TIME column) will
           not be changed.  If the record contains a value >0 in DESEL-
           STRAT,  let  us suppose a 1, then all records with the value
           of 1 in DESELSTRAT will share the same time value (or  what-
           ever  DESEL  selected),   extimated together.  Those records
           with value 2 will be another set of records which will share
           a  time  value, etc.  Thus, within a subject, there may be a
           PK record and a PD record which should share the  same  time
           value.   Also,  a  group of subjects may share the same time
           values.  Within  a  subject,  times  will  be  automatically
           sorted  as they are changed, so that NONMEM's  time ordering
           policy is not violated.

      DESELMIN=label
           The data item (column) containing the minimal value

      DESELMAX=label
           The data item (column) containing the maximal value

           You must impose boundaries on  the  values  that  are  being
           optimized.   That  is  done with these two data items.  Only
           those records with a stratification value >0  in  DESELSTRAT
           column  will  require  a  min  and max value, and only those
           records that define that stratification value for the  first
           time.

           All  four   DESEL  items  must  be  specified:  DESEL,DESEL-
           STRAT,DESELMIN,DESELMAX.  They may be repeated for  as  many
           design elements are to be optimized.
           For example for times and amounts:
           DESEL=TIME DESELMMIN=TMIN  DESELMAX=TMAX   DESELSTRAT=TSTRAT
           DESEL=AMT  DESELMIN=AMTMIN DESELMAX=AMTMAX DESELSTRAT=AMTSTRAT

      NMIN=label
           The name of the data item (column) containing minimal number
           of time points to the subject.

      NMAX=label
           The name of the data item (column) containing maximal number
           of time points to the subject.
           If  NMIN<0  or NMAX<0, then most previous non-negative value
           is used.  The NMIN and NMAX column are  only  used  for  the
           DISCRETE*  analyses, to bound the number of time points that
           may be permitted for a given subject.  With  DISCRETE*,  the
           total  N  of time points among all subjects is determined by
           the total number of time points whose MDV=0 in the  starting
           data set.

      STRAT=label
           The data item (column) containing grouping or stratification
           number pertaining to that subject.

      STRATF=label
           The data item (column) containing starting  fraction  repre-
           sentation for the stratification value in column STRAT.
           If STRAT and STRATF are specified, and there is at least one
           STRAT value >0, then the SRATF  values  are  optimized,  and
           represent  the weight to the contribution of that subject to
           the Information matrix.  For  STRAT<=0,  then  their  STRATF
           values  are not optimized, and remain fixed at their initial
           values, but are still used as  weights  to  the  information
           matrix.   It  is  up  to  the user to ensure that the sum of
           STRATF values among unique STRAT values sum to 1.  If  value
           of  STRATF<0.0,  then  that  subject  is not included in the
           assessment.

 Additional $COVARIANCE Control Options for $DESIGN

 In addition, $DESIGN sets up the  covariance  step  as  $COV  MATRIX=R
 UNCONDITIONAL  without  the  user requiring this record entered in the
 control stream.  If you wish to specify  additional  control  for  the
 covariance, you can add these in a $COV record, such as:
 $DESIGN FIMDIAG=1 OFVTYPE=6 APPROX=FO MAXEVAL=9999 NOHABORT PRINT=20 SIGL=10 POSTHOC
 $COV PRINT=E

 RECORD ORDER:

 Same as $ESTIMATION or $SIMULATION

 REFERENCES: Guide Introduction_7
