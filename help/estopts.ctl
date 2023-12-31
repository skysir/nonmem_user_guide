


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     $ESTIMATION RECORD OPTIONS                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Instructions for the NONMEM Estimation Step
 CONTEXT: NM-TRAN Control Record

 List  of  $ESTIMATION  Record  Options  and Their Relevance to Various
 Methods

 Classical methods are the First order method (METHOD=0 or  METHOD=FO),
 the Conditional method (METHOD=1 or METHOD=FOCE) and the Hybrid method
 (METHOD=HYBRID).  All other methods are EM  (Expectation-Maximization)
 / Bayesian methods and are new to NONMEM 7.

 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS
 -2LL                    X           X     X        X     X        X      X       X
 ATOL +                  X           X     X        X     X        X      X       X
 AUTO                                X     X        X     X        X      X       X
 BAYES_PHI_STORE                                                          X       X
 BIONLY                                                                   X       X
 BOOTDATA                X           X     X        X     X        X
 CALPHA                              X     X        X     X        X      X       X
 CENTERING               X
 CINTERVAL                           X     X        X     X        X      X       X
 CITER/CNSAMP                        X     X        X     X        X      X       X
 CLOCKSEED                                 X        X     X        X      X       X
 CONDITIONAL             X           X     X        X     X        X      X       X
 CONSTRAIN                           X     X        X     X        X      X       X
 CTYPE                   (CTYPE 4)   X     X        X     X        X      X       X
 DERCONT                                   X        X     X        X
 DF                                                 X     X
 DFS(CHAIN only)
 EONLY                                     X        X     X        X
 ETABARCHECK             X
 ETADER                  X           X              X     X
 ETASAMPLES                                                        X      X       X
 ETASTYPE                X           X     X        X     X        X
 FAST                    X           X
 FILE                    X           X     X        X     X        X      X       X
 FNLETA                  X           X     X        X     X        X      X       X
 FORMAT/DELIM            X           X     X        X     X        X      X       X
 FPARAFILE               X           X     X        X     X        X      X       X
 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS
 GRD                                 X     X        X     X        X      X       X
 GRDQ                                               X     X
 GRID(Stieltjes)         X
 HYBRID                  X
 IACCEPT                                            X     X        X      X       X
 IACCEPTL                                           X     X
 INTERACTION             X           X     X        X     X        X      X       X
 IKAPPA                                                                           X
 ISAMPEND                                  X        X     X        X
 ISAMPLE                                                           X      X       X
 ISAMPLE_M1                                                        X      X       X
 ISAMPLE_M1A                                                       X      X       X
 ISAMPLE_M1B                                                       X      X       X
 ISAMPLE_M2                                                        X      X       X
 ISAMPLE_M3                                                        X      X       X
 ISCALE_MAX                                         X     X        X      X       X
 ISCALE_MIN                                         X     X        X      X       X
 KAPPA                                                                            X
 KNUTHSUMOFF             X           X     X        X     X        X      X       X
 LAPLACE                 X           X     *        *     X        *      *       *
 LEVCENTER(for $LEVEL)   X           X     X        X     X        X      X       X
 LEVWT (for $LEVEL)      X           X     X        X     X        X      X       X
 LIKE                    X           X     X        X     X        X      X       X
 LNTWOPI                 X           X     X        X     X        X      X       X
 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS
 MADAPT                                                                           X
 MAPCOV                                             X     X
 MAPINTER                                           X     X
 MAPITER                                            X     X
 MAPITERS                                                          X      X
 MASSRESET               ***                                              X
 MAXEVAL                 X
 MCETA                   X           X              X     X
 MSFO                    X           X     X        X     X        X      X       X
 MUM                                 X     X        X     X        X      X       X
 NBURN                                                             X      X       X
 NITER/NSAMPLE                       X     X        X     X        X      X       X
 NOABORT                 X           X     X        X     X        X      X       X
 NOCOV**                             X     X        X     X        X      X       X
 NOHABORT                X           X     X        X     X        X      X       X
 NOLABEL                 X           X     X        X     X        X      X       X
 NONINFETA               X
 NOOMEGABOUNDTEST        X
 NOPRIOR                 X           X     X        X     X        X      X       X
 NOSIGMABOUNDTEST        X
 NOSLOW                  X           X
 NOSUB                   X           X     X        X     X        X      X       X
 NOTHETABOUNDTEST        X
 NOTITLE                 X           X     X        X     X        X      X       X
 NSIG/SIGDIGITS          X           X
 NUMDER                  X           X     X        X     X        X      X       X
 NUMERICAL               X           X     *        *     X        *      *       *
 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS
 NUTS_BASE                                                                        X
 NUTS_DELTA                                                                       X
 NUTS_EPARAM             ***                                              X
 NUTS_GAMMA                                                                       X
 NUTS_INIT                                                                        X
 NUTS_MASS               ***                                              X
 NUTS_MAXDEPTH                                                                    X
 NUTS_OPARAM             ***                                              X
 NUTS_REG                ***                                              X
 NUTS_SPARAM             ***                                              X
 NUTS_STEPINTER                                                                   X
 NUTS_STEPITER                                                                    X
 NUTS_TERM                                                                        X
 NUTS_TEST                                                                        X
 NUTS_TRANSFORM                                                                   X
 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS
 OACCEPT                                                                  X       X
 OLKJDF                                                                           X
 OLNTWOPI                                                          X      X       X
 OMEGABOUNDTEST          X
 OMITTED                 X           X     X        X     X        X      X       X
 OPTMAP                  X           X              X     X
 ORDER                   X           X     X        X     X        X      X       X
 OSAMPLE_M1                                                               X       X
 OSAMPLE_M2                                                               X       X
 OSAMPLE_M3                                                               X       X
 OVARF                                                                            X
 PACCEPT                                                                  X       X
 PARAFILE                X           X     X        X     X        X      X       X
 PARAFPRINT              X           X     X        X     X        X      X       X
 PHITYPE                             X     X        X     X        X      X       X
 POSTHOC                 X           X     X        X     X        X      X       X
 PREDICTION              X           X     X        X     X        X      X       X
 PRINT                   X           X     X        X     X        X      X       X
 PRIORC                  X           X     X        X     X        X      X       X
 PSAMPLE_M1                                                               X       X
 PSAMPLE_M2                                                               X       X
 PSAMPLE_M3                                                               X       X
 PSCALE_MAX                                                               X       X
 PSCALE_MIN                                                               X       X
 RANMETHOD=nSmP                            X        X     X        X      X       X
 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS
 REPEAT                  X
 REPEAT1                 X
 REPEAT2                 X
 SADDLE_HESS             X
 SADDLE_RESET            X
 SEED                                      X        X     X        X      X       X
 SIGL                    X           X     X        X     X        X
 SIGLO                   X           X              X     X
 SIGMABOUNDTEST          X
 SLKJDF                                                                           X
 SLOW                    X           X
 SORT                    X
 STDOBJ                                             X     X
 STIELTJES               X
 SVARF                                                                            X
 TBLN(CHAIN only)
 THETABOUNDTEST          X
 THIN                                                                     X       X
 TPU                                                                              X
 TTDF                                                                             X
 ZERO                    X
 Option                  Classical   ITS   DIRECT   IMP   IMPMAP   SAEM   BAYES   NUTS

 + ADVAN9,ADVAN13,ADVAN14,ADVAN15,ADVAN16,ADVAN17,ADVAN18
 *May be needed to suppress error messages from NMTRAN or NONMEM.
 **When last estimation step
 ***In prep for NUTS

 (See $ESTIMATION).

 REFERENCES: Guide Introduction_7
