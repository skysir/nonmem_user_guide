


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           PRED,RES,WRES                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied PRED and INFN routines

 USAGE:
      USE ROCM_REAL, ONLY:                                             &
         PRED_=>APPND(001),RES_=>APPND(002),WRES_=>APPND(003)          &
         IWRS_=>APPND(004),IPRD_=>APPND(005),IRS_=>APPND(006)          &
         NPRED_=>APPND(007),NRES_=>APPND(008),NWRES_=>APPND(009)       &
         NIWRES_=>APPND(010),NIPRED_=>APPND(011),NIRES_=>APPND(012)    &
         CPRED_=>APPND(013),CRES_=>APPND(014),CWRES_=>APPND(015)       &
         CIWRES_=>APPND(016),CIPRED_=>APPND(017),CIRES_=>APPND(018)    &
         PREDI_=>APPND(019),RESI_=>APPND(020),WRESI_=>APPND(021)       &
         IWRESI_=>APPND(022),IPREDI_=>APPND(023),IRESI_=>APPND(024)    &
         CPREDI_=>APPND(025),CRESI_=>APPND(026),CWRESI_=>APPND(027)    &
         CIWRESI_=>APPND(028),CIPREDI_=>APPND(029),CIRESI_=>APPND(030) &
         EPRED_=>APPND(031),ERES_=>APPND(032),EWRES_=>APPND(033)       &
         EIWRES_=>APPND(034),EIPRED_=>APPND(035),EIRES_=>APPND(036)    &
         NPDE_=>APPND(037),ECWRES_=>APPND(038),NPD_=>APPND(039)        &
         OBJI_=>APPND(040)

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: APPND(MAXXNAME)

 DISCUSSION:

 The PRED_,RES_,WRES_ items are the same as PRED, RES and WRES items.

 The  CPRED_,CRES_,CWRES_  are calculated as if the estimation was per-
 formed using a CONDITIONAL method without INTERACTION.   (conditional,
 non-interactive in the manner of Hooker et al.)

 The  PREDI_,RESI_,WRESI_  are calculated as if the estimation was per-
 formed using the FO method with INTERACTION (non-conditional, interac-
 tive).

 The  CPREDI_,CRESI_,CWRESI_  are  calculated  as if the estimation was
 performed using a CONDITIONAL method  with  INTERACTION  (conditional,
 interactive)

 EPRED_,ERES_,EWRES_ are Monte-Carlo generated pred, res, and wres val-
 ues, and are not linearized approximations  like  the  other  weighted
 residual types.

 NPDE_  is  the  normalized  prediction  distribution error (takes into
 account within-subject correlations)

 NPD_ is the correlated normalized prediction distribution error  (does |
 not take into account within-subject correlations)

 EWRES_  and NPDE_  and NPD_ are evaluated using predicted function and
 residual variances evaluated over a Monte Carlo sampled range of  etas
 with population variance Omega.

 ECWRES_ is evaluated with only the predicted function evaluated over a
 Monte Carlo sampled range of  etas  with  population  variance  Omega,
 while  residual  variance  is  always evaluated at conditional mode or
 mean.  Thus, ECWRES is the Monte Carlo version of CWRES,  while  EWRES
 is the Monte Carlo version of CWRESI.

 The  EPRED_,  ERES_, EWRES_, NPD_, NPDE_, ECWRES_ items are calculated |
 by the Table Step.  When EPRED_, ERES_, EWRES_, NPD_,  NPDE_,  ECWRES_ |
 items are displayed, the options ESAMPLE and SEED of the $TABLE record |
 are of interest.  (See $table).

 The PRED_,RES_,WRES_ (PRED,RES,WRES) items behave  as  they  did  with
 NONMEM  VI. Consequently, if FO or a CONDITIONAL method without INTER-
 ACTION is used, NPRED_,NRES_,NWRES_ and PRED_,RES_,WRES_  are  equiva-
 lent.   If  FO  or  a  CONDITIONAL  method  with  INTERACTION is used,
 IPRED_,IRES_,IWRES and PRED_,RES_,WRES_ are equivalent.

 These items are available with passes  through  the  data  set  during
 problem finalization (i.e. ICALL=3), in a data record specific way.

 These  items  may  be used as right-hand quantities in $PRED and $INFN
 during problem finalization.  They can also be output in  table  files
 by  specifying their names (without the underscore) on a $TABLE state-
 ment.

 (See pred res, wres).

 Location prior to NONMEM 7: rocm49

 For other variables, see the INTRODUCTION TO NONMEM 7.

 REFERENCES: Guide I, section C.3.5.3 , C.3.5.4 
 REFERENCES: Guide IV, section III.B.16 , III.B.17 
 REFERENCES: Guide V, section 9.5 , 10.7 , 11.4.4.2  
 REFERENCES: Guide Introduction_7
