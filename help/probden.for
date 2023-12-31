


 +--------------------------------------------------------------------+
 |                                                                    |
 |                   PROBABILITY DENSITY FUNCTIONS                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Functions that may be used in abbreviated code.
 CONTEXT: Fortran coded function

 A  series of built in probability density functions are available with
 NONMEM 7.4.2.  For a given probability density there is also a cumula-
 tive  distribution function (densitycdf), and random number generating
 function (density_rng).

 They are described in INTRODUCTION TO NONMEM 7.4.2 Section I.26.

 The density and densityCDF functions have arguments that are  compati-
 ble with the FUNC system, in which function provides derivatives (XD),
 and second derivatives (XDD) (see I.65.Expanded  Syntax  and  Capacity
 for  User-Defined  Functions  (FUNCA)  (NM74)). Thus, even random (eta
 associated) variables may serve as arguments to the parameters of  the
 density   functions.  The  source  code  of  these  densities  are  in
 ..\source\DISTRIB.f90, DISTRIBCDF.f90, and DISTRIBRNG.f90.  Note  that
 multi-variate densities do not have a corresponding CDF routine.

 Examples of their use are in the ..\examples\densities directory.

 Here  are  the  list  of densities, which are modeled after the format
 from the Stan manual:

 BERNOULLI
 BERNOULLILOGIT
 BINOMIAL
 BINOMIALLOGIT
 BETABINOMIAL
 HYPERGEOMETRIC
 CATEGORICAL
 CATEGORICALLOGIT
 ORDEREDLOGISTIC
 NEGBINOMIAL
 NEGBINOMIAL2
 NEGBINOMIAL2LOG
 POISSON
 POISSONLOG
 MULTINOMIAL
 NORMAL
 EXPMODNORMAL
 SKEWNORMAL
 STUDENTT
 DOUBLEEXPONENTIAL
 LOGISTIC
 GUMBEL
 LOGNORMAL
 CHISQUARE
 INVCHISQUARE
 SCALEDINVCHISQUARE
 EXPONENTIAL
 GAMMA
 INVGAMMA
 WEIBULL
 FRECHET
 RAYLEIGH
 PARETO
 PARETO2
 BETA
 DIRICHLET
 VON MISES

 REFERENCES: none
