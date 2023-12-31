


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              MU_MODEL                              |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING:
 CONTEXT: Abbreviated code

 USAGE:
 MU_3=LOG(THETA(3))
 V=EXP(MU_3+ETA(3))

 DISCUSSION:

 The new methods in NONMEM are most efficiently implemented if the user
 supplies information on how the THETA parameters are associated arith-
 metically  with  the  etas  and individual parameters, wherever such a
 relationship holds.  Calling the individual parameters phi, the  rela-
 tionship should be
 phi_i=mu_i(theta)+eta(i)
 for each parameter i that has an eta associated with it, and mu_i is a
 function of THETA.  The association of one or more THETA's with ETA(1)
 must be identified by a variable called MU_1.  Similarly, the associa-
 tion with ETA(2) is MU_2, that of ETA(5) is MU_5, etcetera.   This  is
 called "MU Referencing", or "MU Modelling".

 Providing  this  information is as straight-forward as introducing the
 MU_ variables into the $PRED or $PK code by expansion of the code.

 For  a very simple example, the original code may have the lines

 CL=THETA(4)+ETA(2)

 This may be rephrased as:

 MU_2=THETA(4)
 CL=MU_2+ETA(2)

 Another example would be:

 CL=(THETA(1)*AGE**THETA(2))*EXP(ETA(5))
 V=THETA(3)*EXP(ETA(3))

 which would now be broken down into two  additional  lines,  inserting
 the definition of a MU as follows:

 MU_5= LOG(THETA(1))+THETA(2)*LOG(AGE)
 MU_3=LOG(THETA(3))
 CL=EXP(MU_5+ETA(5))
 V=EXP(MU_3+ETA(3))

 Note  the  arithmetic  relationship  identified by the last two lines,
 where MU_5+ETA(5) and MU_3+ETA(3) are expressed. This action does  not
 change the model in any way.

 If  the  model  is  formulated  by  the traditional typical value (TV,
 mean), followed by individual value, then it  is  straight-forward  to
 add the MU_ references as follows:

 TVCL= THETA(1)*AGE**THETA(2)
 CL=TVCL*EXP(ETA(5))
 TVV=THETA(3)
 V=TVV*EXP(ETA(3)
 MU_3=LOG(TVV)
 MU_5=LOG(TVCL)

 This  also  will work because only the MU_x= equations are required in
 order to take advantage of EM efficiency.  It is not required  to  use
 the  MU_  variables in the expression EXP(MU_5+ETA(5)), since the fol-
 lowing are equivalent:
 CL=TVCL*EXP(ETA(5))=EXP(LOG(TVCL)+ETA(5)=EXP(MU_5+ETA(5))
 but it helps as an exercise to determine that the  MU_  reference  was
 properly  transformed (in this case log transformed) so that it repre-
 sents an arithmetic association with the eta.

 An incorrect usage of MU modeling would be:

 MU_1=LOG(THETA(1))
 MU_2=LOG(THETA(2))
 MU_3=LOG(THETA(3))
 CL=EXP(MU_1+ETA(2))
 V=EXP(MU_2+MU_3+ETA(1))

 In the above example, MU_1 is used as an arithmetic  mean  to  ETA(2),
 and  a  composite  MU_2  and  MU_3 are the arithmetic means to ETA(1),
 which would not be correct.  The association of  MU_x+ETA(x)  must  be
 strictly adhered to.

 Once one or more thetas are modeled to a MU, the theta may not show up
 in any subsequent lines of code.  That is,  the  only  usage  of  that
 theta may be in its connection with MU.  For example, if

 CL=THETA(5)*EXP(ETA(2))

 it can be rephrased as

 MU_2=LOG(THETA(5))
 CL=EXP(MU_2+ETA(2))

 But  later,  suppose  THETA(5)  is  used  without its association with
 ETA(2):

  ...
 CLZ=THETA(5)*2

 Then THETA(5) cannot be MU modeled, because it shows up as  associated
 with  ETA(2) in one context, but as a fixed effect without association
 with ETA(2) elsewhere.

 However, if

 MU_2=LOG(THETA(5))
 CL=EXP(MU_2+ETA(2))
 CLZ=CL*2

 Then this is legitimate, as the individual parameter  CL  retains  the
 association  of  THETA(5)  with ETA(2), when used to define CLZ.  That
 is, THETA(5) and ETA(2) may not used separately in any other  part  of
 the model, except indirectly through CL, in which their association is
 retained.

 Suppose you have:

 CL=THETA(5)+THETA(5)*ETA(2)

 One should see this as:

 CL=THETA(5)*(1+ETA(2))

 So the way to MU model this is:

 MU_2=1.0
 CL=THETA(5)*(MU_2+ETA(2))

 Which would mean that in the end, THETA(5) is not actually MU modeled,
 since MU_2 does not depend on THETA(5).  One would be tempted to model
 as follows:

 MU_2=THETA(5)
 CL=MU_2+MU_2*ETA(2)

 But this would be incorrect, as  MU_2  and  ETA(2)  may  not  show  up
 together  in  the code except as MU_2+ETA(2) or its equivalent.  Thus,
 THETA(5) cannot be MU modeled.  In such cases, remodel to the  follow-
 ing similar format:

 CL=THETA(5)*EXP(ETA(2))

 So that THETA(5) may be MU modeled as:

 MU_2=LOG(THETA(5))
 CL=EXP(MU_2+ETA(2))

 Sometimes,  a  particular  parameter has a fixed effect with no random
 effect, such as:

 Km=THETA(6)

 with the intention that Km is unknown but  constant  across  all  sub-
 jects.   In  such  cases, the THETA(6) and Km cannot be Mu referenced,
 and the EM efficiency will not be  available  in  moving  this  Theta.
 However,  one  could assign an ETA to THETA(5), and then fix its OMEGA
 to a small value, such as 0.0225 =0.15^2 to represent 15% CV, if OMEGA
 represents  proportional  error.   This  often will allow the EM algo-
 rithms to efficiently move this parameter, while retaining the  origi-
 nal  intent  that  all  subjects have similar, although not identical,
 Km's.  Very often, inter-subject variances to parameters were  removed
 because  the  FOCE had difficulty estimating a large parametered prob-
 lem, and so it was an artificial constraint to begin with.  EM methods
 are  much  more  robust,  and  are adept at handling large, full block
 OMEGA's, so you may want to incorporate as many etas as possible  when
 using the EM methods.

 You  should  Mu  reference  as many of the THETA's as possible, except
 those pertaining to residual variance (which should be modeled through
 SIGMA  whenever  possible).   If you can afford to slightly change the
 theta/eta relationship a little  to  make  it  MU  referenced  without
 unduly  influencing the model specification or the physiological mean-
 ing, then it should be done.

 When the arithmetic mean of an ETA is  associated  with  one  or  more
 THETA's in this way, EM methods can more efficiently analyze the prob-
 lem, by requiring in certain calculations only the evaluation  of  the
 MU's  to  determine  new  estimates  of THETAs for the next iteration,
 without having to re-evaluate the predicted value  for  each  observa-
 tion,  which  can be computationally expensive, particularly when dif-
 ferential equations are used in the model.  For those THETA's that  do
 not  have  a  relationship  with any ETA's, and therefore cannot be MU
 referenced (including THETA's associated with ETAS whose  OMEGA  value
 is fixed to 0), computationally expensive gradient evaluations must be
 made to provide new estimates of them for the next iteration.

 There is additional increased efficiency  in  the  evaluation  of  the
 problem  if  the MU models are linear functions with respect to THETA.
 Recalling one of the previous examples above, we could re-parameterize
 THETA such that

 MU_5=THETA(1)+THETA(2)*LOG(AGE)
 CL=EXP(MU_5+ETA(5))
 MU_3=THETA(3)
 V=EXP(MU_3+ETA(3))

 This  changes  the  values  of THETA(1) and THETA(3) such that the re-
 parameterized THETA(1) and THETA(3) are the logarithm of the  original
 parameterization  of THETA(1) and THETA(3).  The models are identical,
 however, in that the same maximum likelihood value will  be  achieved.
 The  only  inconvenience  is  having  to anti-log these THETA's during
 post-processing.

 The added efficiency  obtained  by  maintaining  linear  relationships
 between  the  MU's  and THETA's is greatest when using the SAEM method
 and the MCMC Bayesian method.  In the Bayesian  method,  THETA's  that
 are  linearly  modeled with the MU variables have linear relationships
 with respect to the inter-subject variability,  and  this  allows  the
 Gibbs  sampling  method  to be used, which is much more efficient than
 the Metropolis-Hastings (M-H) method.  By default,  NONMEM  tests  MU-
 THETA  linearity  by  determining  if the second derivative of MU with
 respect to THETA is nearly or equal to 0.  Those THETA parameters with
 0  valued second derivatives are Gibbs sampled, while all other THETAS
 are M-H sampled.  In the Gibbs sampling method, THETA values are  sam-
 pled  from a multi-variate normal conditional density given the latest
 PHI=MU+ETA values  for  each  subject,  and  the  samples  are  always
 accepted.   In  M-H  sampling,  the  sampling  density used is only an
 approximation, so the sampled THETA values must be tested by  evaluat-
 ing the

 Some additional rules for MU referencing are as follows:

 1)   As  much  as  possible, define the MU's in the first few lines of
      $PK or $PRED.  Do not define MU_ values in $ERROR.  Have all  the
      MU's  particularly  defined  before any additional verbatim code,
      such as write statements.  NMTRAN produces a MUMODEL2  subroutine
      based  on  the  PRED  or  PK  subroutine  in  FSUBS.F90, and this
      MUMODEL2 subroutine is frequently called with  the  ICALL=2  set-
      tings,  more  often  than  PRED or PK.  The fewer code lines that
      MUMODEL2 has to go through to evaluate all  the  MU_s'  the  more
      efficient.

 2)   Whenever possible, have the MU variables defined unconditionally,
      outside IF...THEN blocks.

 3)   Time dependent covariates cannot be part  of  the  MU_  equation.
      For example

      MU_3=THETA(1)*TIME+THETA(2)

      should not be done.  Or, consider

      MU_3=THETA(2)/WT

      Where  WT  varies  with  time.   This would also not be suitable.
      However, we could phrase as

      MU_3=THETA(2)
      CL=WT*(MU_3+ETA(3))

      is fine, where MU_3 represents a population  mean  clearance  per
      unit  weight,  which  is  constant  with time, and more universal
      among subjects, whereas CL is the  non-wieght  normalized  clear-
      ance,  than  depends  on a person's weight, which could vary with
      time as well.  The MU variables may vary with inter-occasion, but
      not with time.

 4)   Starting with NONMEM 7.2, NMTRAN's CHECKMU subroutine attempts to |
      look for errors in MU modeling.  If it appears that there may  be |
      errors, then there are messages such as                           |

      (MU_WARNING 13) MU_001: DOES NOT HAVE ADDITIVE ASSOCIATION WITH ETA(001)|

      Such  warnings  do  not  affect the outputs from NMTRAN. FSUBS is |
      generated as usual.  Sometimes the warnings may be  ignored  (see |
      "Model parameters as log t-Distributed", below.)  Sometimes warn- |
      ings may not be generated when they should be.   Thus,  the  user |
      must pay close attention to following the rules.                  |

      Option  NOCHECKMU  of the $ABBR record may be used to prevent NM- |
      TRAN from attempting to check the MU model statements.

 Examples show examples of  MU  modeling  for  various  problem  types.
 Study  these examples carefully. When transposing your own code, begin
 with simple problems and work your way to more complex problems.

 At this point one may wonder why bother  inserting  MU  references  in
 your  code.  MU referencing only needs to be done if you are using one
 of the new EM or Gibbs sampling methods to improve  their  efficiency.
 The  EM methods may be performed without MU references, but it will be
 several fold slower than the FOCE method, and the problem may not even
 optimize  successfully.  For simple two compartment models, the new EM
 methods are slower than FOCE even with the MU references.  But, for  3
 compartment models, or numerical integration problems, the improvement
 in speed by the EM methods, properly MU  modeled,  can  be  5-10  fold
 faster than with FOCE.

 Example  6  described  at  the  end of the SIGL section is one example
 where importance sampling solves this problem in 30  minutes,  with  R
 matrix  standard  error, versus FOCE which takes 2-10 hours or longer,
 and without even requesting the $COV  step.   So,  for  complex  PK/PD
 problems  that take a very long time in FOCE, it is well worth putting
 in MU references and using one of the EM methods, even if you may need
 to  rephrase  some  of  the fixed/random (theta/eta) effects relation-
 ships.  In addition, FOCE is a linearized optimization method, and  is
 less accurate than the EM and Bayesian methods when data are sparse or
 when the posterior density for each individual is highly non-normal.

 Model parameters as log t-Distributed in the Population                |

 Sometimes one may suspect that PK/PD model parameters are actually log |
 t-distributed  among  the  population,  with  degrees  of  freedom NU, |
 instead of the usual log normal distributed.                           |

 See INTRODUCTION TO NONMEM 7, MU Referencing                           |

 An example of simulation and analysis of such data is given as         |

  ..\examples\tdist6_sim.ctl                                            |
  ..\examples\tdist6.ctl:                                               |
  ..\examples\tdist7.ctl                                                |

 Note that constructions such as                                        |
 CL=EXP(MU_1+ETA(1)*SQRT((EXP(CLR)-1.0)/CLR))                           |
 violate the strict  MU_x+ETA(x)  rule  recommended  for  EM  analysis, |
 because  the  term  SQRT((EXP(CLR)-1.0)/CLR)  is multiplied by ETA(1). |
 NM-TRAN will generate a number of  MU_WARNING  messages.   Nonetheless |
 for  this  example,  the importance sampling works quite well, and the |
 MU_WARNING messages may be ignored.

 REFERENCES: Guide Introduction_7
