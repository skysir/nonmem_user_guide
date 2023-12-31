


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                TOL                                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: TOL subroutine
 CONTEXT: User-supplied subroutine; for use with PREDPP

 USAGE:

 Versions before NONMEM 7.4:

 SUBROUTINE TOL(NRD)
 USE SIZES,     ONLY: ISIZE
 INTEGER(KIND=ISIZE) :: NRD
 DIMENSION :: NRD(*)

 With NONMEM 7.4:

 SUBROUTINE TOL(NRD,ANRD,NRDC,ANRDC)
 USE SIZES,     ONLY: ISIZE
 INTEGER(KIND=ISIZE) :: NRD(0:*), ANRD(0:*), NRDC(0:*), ANRDC(0:*)

 Optional declarations with NONMEM 7.4:

 USE NMPRD_INT, ONLY: IPROB
 USE NM_BAYES_INT, ONLY: NM_STEP,BASE_STEP,EST_STEP,COV_STEP, &
  TABLE_STEP,SIML_STEP,INE_STEP, NONP_STEP

 DISCUSSION:
 The TOL subroutine is called by PREDPP when ADVAN 6, 8, 9, 10, 13, 14,
 15, 16, 17, or 18 is used.  It is also called when SS6 or SS9 is used.
 With  NONMEM  7.4,  there  are multiple calls during the run, for each
 NONMEM step.  With earlier version, TOL is called  only  once  at  the
 start of a run.

 Output argument:

  NRD(I)
      The number of digits that are required to be accurate in the com-
      putation of the drug amount in compartment I, i.e., the  relative
      tolerance.   ADVAN  9,  13, 16, 17, and 18 have the capability of
      using
       different values of NRD for different compartments (ADVAN 14 and
      15  only  allow  different  absolute tolerances for each compart-
      ment).  For compartments not specified, the tolerance of the last
      compartment specified will be used.

      However, all the other ADVAN routines requiring TOL take the rel-
      ative tolerance to be the same for all compartments; NRD(I), I  >
      1,  is ignored, and only NRD(1) is used.  With NONMEM 7.4, NRD(0)
      is the relative tolerance for the Steady State computations.   If
      NRD(0) is not specified, NRD(1) is used.

      The  value  of  NRD(1)  can  also be specified using $SUBROUTINES
      option TOL.
      The value of NRD(0) can  also  be  specified  using  $SUBROUTINES
      option SSTOL.

  ANRD(I) (NM74)
      The  absolute  tolerance in the computation of the drug amount in
      compartment  I.   The  default  is  12  (that  is,  accuracy   is
      10**(-12)).   Used  by ADVAN 9, 13, 14, 15, 16, 17, and 18, which
      have the capability of using
       different values of ANRD for different compartments.   For  com-
      partments  not  specified,  the tolerance of the last compartment
      specified will be used.

      ANRD(0) is the absolute tolerance for the Steady  State  computa-
      tions.  If ANRD(0) is not specified, ANRD(1) is used.

      The  value  of  ANRD(1)  can also be specified using $SUBROUTINES
      option ATOL.
      The value of ANRD(0) can also  be  specified  using  $SUBROUTINES
      option SSATOL.

  NRDC(I) (NM74)
      Same  as  NRD(I),  but used for the FOCE/LAPLACE covariance step.
      Used with ADVAN 9, 13, 14, and 15, 16, 17, and 18.  If  not  set,
      NRDC  defaults  to  the value of NRD.  NRDC(0) is used for Steady
      State computations during the FOCE/LAPLACE covariance step.

      The value of NRDC(1) can also  be  specified  using  $SUBROUTINES
      option TOLC.
      The  value  of  NRDC(0)  can also be specified using $SUBROUTINES
      option SSTOLC.

  ANRDC(I) (NM74)
      Same as ANRD(I), but used for the FOCE/LAPLACE  covariance  step.
      Used  with  ADVAN 9, 13, 14, and 15, 16, 17, and 18.  If not set,
      ANRDC defaults to the value of ANRD.  ANRDC(0) is used for Steady
      State computations during the FOCE/LAPLACE covariance step.

      The  value  of  ANRDC(1) can also be specified using $SUBROUTINES
      option ATOLC.
      The value of ANRDC(0) can also be  specified  using  $SUBROUTINES
      option SSATOLC.

 When  NM-TRAN  is used, this information may be supplied by either the
 TOL option of the $SUBROUTINES record or the $TOL record.

 Finally, you may supply a TOL routine that assigns values of  NRD  and
 ANRD  specifically for the initial (base) setting and each NONMEM step
 (estimation, covariance, simulation, table/scatter  step,  simulation,
 initial  parameters  estimate,  nonparametric).  For example, create a
 toluser.f90 file,

 SUBROUTINE TOL(NRD,ANRD,NRDC,ANRDC)
 USE NMPRD_INT, ONLY: IPROB
 USE NM_BAYES_INT, ONLY: NM_STEP,BASE_STEP,EST_STEP,COV_STEP, &
  TABLE_STEP,SIML_STEP,INE_STEP, NONP_STEP
 IMPLICIT NONE
 INTEGER :: NRD(0:*), ANRD(0:*), NRDC(0:*), ANRDC(0:*)
 IF(NM_STEP==EST_STEP) THEN
 NRD(1)=6
 ANRD(1)=10
 ELSE IF (NM_STEP==COV_STEP) THEN
 NRD(1)=7
 ANRD(1)=8
 ELSE IF (NM_STEP==TABLE_STEP) THEN
 NRD(1)=8
 ANRD(1)=7
 ELSE
 NRD(1)=9
 ANRD(1)=12
 ENDIF
 IF(IPROB>1) THEN
 NRD(1)=NRD(1)+1
 ANRD(1)=ANRD(1)+1
 ENDIF
 RETURN
 END

 and incorporate using $SUBR, e.g.,

 $SUBROUTINES ADVAN13 TRANS1 TOL=toluser.f90

 REFERENCES: Guide IV, section V.C.4 
 REFERENCES: Guide VI, section VI.D , Figure 41
 REFERENCES: Guide VI, Appendix II 
