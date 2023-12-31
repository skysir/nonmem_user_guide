


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           NONMEM MODULES                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Global variables in NONMEM
 CONTEXT: User-supplied routines

 DISCUSSION:
 FORTRAN  modules are used to communicate values between various compo-
 nents of the NONMEM system.  They supplement the subroutine arguments.

 1) NONMEM-PRED modules

 These modules contain values that are (for the most part) communicated
 from  PRED  to  NONMEM.   Prior  to  NONMEM 7 these values were COMMON
 blocks.  Following is a list of the MODULE, the (old)  COMMON,  and  a
 description of the variables.

 NMPRD_INT (NMPRD1)
      PRED return code and user message count (IERPRD,NETEXT)

 NMPRD_CHAR (NMPRD2)
      PRED user message (ETEXT)

 NMPRD_INT (NMPRD3)
      "Copying  pass" flag and SAVE region size for NMPRD4 (COMACT,COM-
      SAV)

 NMPRD4P (NMPRD4)
      PRED-defined items for tables and scatterplots

 NMPRD_REAL (NMPRD5)
      Correlation matrix for epsilon's.

 NMPRD_INT (NMPRD6)
      Return code from routine SIMEPS.

 NMPRD_REAL (NMPRD7)
      Simulated etas for tables and scatterplots.

 NMPRD_REAL (NMPRD8)
      PRED's "recursive" flag

 NMPRD_REAL (NMPRD9)
      Data record at ICALL 0, 1 and 3

 NMPR_INT (NMPR10)
      Control information for PRED repetition feature.

 NMPR_INT (NMPR11)
      Flags to override effect of NEW option on $SIMULATION record.

 NMPR_REAL (NMPR12)
      Conditional limits on observation.

 NMPR_REAL (NMPR13)
      Lower limits defining interval datum, and their derivatives.

 NMPR_REAL (NMPR14)
      Upper limits defining interval datum, and their derivatives.

 NMPR_INT (NMPR15)
      Skip control variable.

 NMPR_REAL (NMPR16)
      Parameter values produced during Simulation Step.

 NMPR_INT (NMPR17)
      Flag indicating character of PRED variable F.

 2) NONMEM read-only modules

 These modules contain read-only values that are communicated from NON-
 MEM to PRED and other user-supplied routines.

 ROCM_REAL (ROCM0)
      Current theta (THETA)

 ROCM_REAL (ROCM1)
      DV, data items in current L1 record

 ROCM_INT (ROCM1C)
      Size  of  current L1 record (prior to NONMEN 7, this was combined
      with ROCM1).

 ROCM_INT (ROCM2)
      Number of L2 records and length of L2 record

 ROCM_REAL (ROCM3)
      Predictions and derivatives with current L1 record

 ROCM_REAL (ROCM4)
      Selected data from an individual record

 ROCM_INT (ROCM4)
      Selected data from an individual record

 ROCM_REAL (ROCM5)
      Predictions and derivatives with current L2 record

 ROCM_REAL (ROCM6)
      Initial/final theta, omega, sigma

 ROCM_REAL (ROCM7)
      Standard errors of estimates of theta, omega, sigma

 ROCM_REAL (ROCM8)
      Final value of the objective function

 ROCM_INT (ROCM9)
      Return codes from Estimation and Covariance Steps

 ROCM_INT (ROCM10)
      Simulation: no. of replications (total and current)

 ROCM_INT (ROCM11)
      Mixture: index of subpopulation (current and most probable)

 ROCM_INT (ROCM12)
      Second and first derivative flags for PRED

 ROCM_INT (ROCM13)
      Simulation: Final pass flag for PRED and CONTR

 NMPRD_INT (ROCM14)
      Problem and subproblem counters.

 NMPRD_INT (ROCM15)
      List of inactive etas.

 ROCM_REAL (ROCM16)
      The number of significant digits in the the vector of final esti-
      mates.

 ROCM_INT (ROCM17)
      New level 2 record flag

 ROCM_REAL (ROCM18)
      Nonparametric   estimates:   height  and  heights  of  cumulative
      marginals

 ROCM_REAL (ROCM22)
      Current value of omega.

 ROCM_REAL (ROCM25)
      Mixing probabilities

 ROCM_REAL (ROCM28)
      Superproblem printing indicator

 NMPRD_INT (ROCM29)
      Population vs. single-subject data flag.

 ROCM_REAL (ROCM31)
      Template record at ICALL=6.

 ROCM_INT (ROCM32)
      Number of the individual record and of the data record.

 ROCM_INT (ROCM34)
      NEWIND at ICALL 0, 1 and 3.

 ROCM_INT (ROCM35)
      Numbers of thetas, etas and epsilons in the problem.

 ROCM_REAL (ROCM36)
      Individual's  posterior  variance-covariance  matrix.   No   help
      entry.

 ROCM_INT (ROCM37)
      Indicator that Estimation Step is omitted and obj. func. is being
      evaluated  for first time.  No help entry.

 ROCM_REAL (ROCM38)
      Conditional limits for observations  in  individual  record.   No
      help entry.

 ROCM_REAL (ROCM39)
      Conditional limits for observation.  No help entry.

 ROCM_REAL (ROCM40)
      Lower  Limits  defining  interval datum, and their derivatives in
      individual record. No help entry.

 ROCM_REAL (ROCM41)
      Upper Limits defining interval datum, and  their  derivatives  in
      individual record. No help entry.

 ROCM_REAL (ROCM42)
      Lower  Limits defining interval datum, and their derivatives.  No
      help entry.

 ROCM_REAL (ROCM43)
      Upper Limits defining interval datum, and their derivatives.   No
      help entry.

 ROCM_REAL (ROCM44)
      Probability that observation is within (outside) limits.

 ROCM_REAL (ROCM45)
      Probability that category occurs.

 ROCM_INT (ROCM46)
      Number  of  individual  records with observations, and indices of
      first and last such record.

 ROCM_REAL (ROCM47)
      Mixture probabilities with individual record containing  template
      record.

 ROCM_INT (ROCM48)
      Length of individual record.

 ROCM_REAL (ROCM49)
      Prediction,  residual, and weighted residual values for all meth-
      ods, as well as normalized probability distribution error.

 ROCM_INT (ROCM50)
      ID data item for the individual contributions  to  the  objective
      function.

 ROCM_REAL (ROCM50)
      Individual contributions to the objective function.

 Except  as  noted,  each has its own entry in the Help document.  Some
 modules are not fully documented.  The interested user may be able  to
 obtain more information by studying the appropriate sections of NONMEM
 code and previous examples that may be available from advanced users.

 (See PREDPP_modules).

 REFERENCES: Guide IV, section IV.E 
 REFERENCES: Guide VI, section III.E.4 , III.I ,  Figures 5, 6
