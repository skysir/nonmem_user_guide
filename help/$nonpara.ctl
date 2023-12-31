


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           $NONPARAMETRIC                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Instructions for the NONMEM Nonparametric Step
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $NONPARAMETRIC [MARGINALS|ETAS] [MSFO=filename] [RECOMPUTE]
                [EXPAND] [NPSUPP=n | NPSUPPE=n]
                [BOOTSTRAP [STRAT=label] [STRATF=label]]
                [PARAFILE=[filename|ON|OFF]                             |
                [NPESTIM=[0|1]] [NPMAXITER=n]                           |
                [UNCONDITIONAL|CONDITIONAL]  [OMITTED]

 SAMPLE:
 $NONPARAMETRIC    ETAS

 DISCUSSION:
 Optional.  Requests that the NONMEM Nonparametric Step be implemented.
 When present, the $ESTIMATION record must also  be  present  and  must
 specify METHOD=1 or POSTHOC.

 For a given eta, the points of support are the vector of posthoc esti-
 mates of that eta for all individuals (i.e., the CPE values  for  that
 eta),  which  is essentially equal to the number of individuals in the
 data set.

 OPTIONS:

 MARGINALS
      Requests that marginal cumulatives  be  obtained  (the  default).
      These values are found in NONMEM global variables
      (See Nonparametric Density: DEN_,CDEN_)

 ETAS Requests that conditional (nonparametric) estimates of eta values
      be obtained. (Also called the CNPE).

 MSFO=filename
      A Model Specification File  is  output  to  the  given  filename.
      Filename  may  not contain embedded spaces.  If filename contains
      commas, semicolons, equal  signs,  or  parentheses,  it  must  be
      enclosed  in  quotes  ('  or ").  Filename may contain at most 71
      characters.  If filename is the same as any option of  the  $NON-
      PARAMETRIC  record,  it  must  be  enclosed  in quotes.  The MSFO
      option may appear without a file name.  More  precisely,  if  the
      $ESTIMATION  record  is  also  present, and it also specifies the
      MSFO option, then the filename is required only on one of the two
      records,  $ESTIMATION  or  $NONPARAMETRIC,  whichever one appears
      first in the control stream.  If the filename appears present  on
      both  records, it must be the same on both records.  If the file-
      name is omitted on the second of the two records, the MSF  option
      must be the final option on that record.  Default: No MSF is out-
      put.

 RECOMPUTE
      Requests that the nonparametric density estimate occurring in  an
      input MSF should be ignored; the nonparametric estimate should be
      recomputed.

 EXPAND
      After the parametric estimation is performed, the final  eta  MAP |
      (or empirical Bayes estimates, EBE) estimates, based on the final |
      SIGMAS, OMEGAS, and THETAS, are normally used as support  points. |
      When EXPAND is selected, an alternative set of EBEs are evaluated |
      using the initial OMEGA values, but using the  final  THETAS  and |
      SIGMAS.                                                           |

 NPSUPP=n                                                               |
      Number  of  total support points to be used.  If NPSUPP>number of |
      subjects, then extra support points are randomly created from the |
      final  OMEGAS (even when EXPAND is selected).  Only one of NPSUPP |
      or NPSUPPE may be specified.                                      |

 NPSUPPE=n                                                              |
      Number of total support points to be used.  If NPSUPPE>number  of |
      subjects, then extra support points are randomly created from the |
      initial, presumably inflated, OMEGAS (even  when  EXPAND  is  not |
      selected).  Only one of NPSUPP or NPSUPPE may be specified.       |

 BOOTSTRAP                                                              |
      The  original data set is fitted during the parametric estimation |
      ($EST), and the eta support points from the original data set are |
      used for the nonparametric version.  However, a bootstrap sample, |
      with subjects uniformly randomly selected with  replacement  from |
      the original data set, is used for the nonparametric distribution |
      analysis.                                                         |

 STRAT                                                                  |
      The label of a data item that serves as the stratification.  This |
      splits  the  data set into distinct sub-sets, guaranteeing a spe- |
      cific number of subjects will be selected from each category.     |

 STRATF                                                                 |
      The label of a data item that contains the fraction  that  should |
      represent  a  category  in  the  bootstrapped  data set.  Without |
      STRATF, the number of subjects to be taken from a given  category |
      is proportional to the number of subjects.

 CONDITIONAL
      The  Nonparametric  Step  is implemented only when the Estimation
      Step terminates successfully or is  not  implemented  (i.e.,  the
      $ESTIMATION record specifies MAXEVAL=0).  This is the default.

 NPESTIM=[0|1]
      The  default  non-parametric estimation method for assessing sup-
      port point probabilities is an (non-Monte Carlo) expectation-max-
      imization  (EM)  method.   You  may  also choose the non-negative
      least squares method (NNL, [29]), with NPESTIM=1.   While  it  is
      touted  to  be  faster  than  EM (NNL is quadratically convergent
      whereas EM is linearly convergent), several tests have not  indi-
      cated  that  NNL  is  any strong speed advantage.  The reason is,
      that the computation time increases by at least the square of the
      number of support points (MAX(NSPSUPP,NIND)) with the NNL method,
      (least squares methods require  matrix  inversion,  which  is  at
      least  an N2 order process), whereas with EM the computation time
      increases in proporation to MAX(NSPSUPP,NIND).  Thus  the  larger
      the  number of support points, the greater the speed advantage of
      the EM method.

 NPMAXITER=n
      The default maximum iterations for non-parametric  estimation  of
      assessing  support  point probabilities is 1000, which is usually
      more than enough.

 PARAFILE=filename
      As of NONMEM 7.4, nonparametric analysis can be parallelized.     |
      PARAFILE=filename specifies a different parafile  than  was  used |
      for the previous step.                                            |
      PARAFILE=ON turns on parallelization for the Nonparametric Step.  |
      PARAFILE=OFF  turns  off  parallelization  for  the Nonparametric |
      Step.

 UNCONDITIONAL
      The Nonparametric Step is always implemented.

 OMITTED
      The Nonparametric Step is not implemented.

 REFERENCES: None.
