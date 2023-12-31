


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                $TOL                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Marks the beginning of abbreviated code for the TOL routine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $TOL
 abbreviated code

 DISCUSSION:
 The  $TOL  record  is used to specify compartment-specific NRD values.
 It is used with PREDPP's general non-linear  models  (ADVAN6,  ADVAN8,
 ADVAN9,  ADVAN13, ADVAN14, ADVAN15, ADVAN16, ADVAN17, ADVAN18, and SS6
 and SS9).  NRD stands for "Number of Required  Digits,"  although  the
 precise  meaning  depends  on  the particular ADVAN or SS routine that
 uses it.

 $TOL cannot be present if the $SUBROUTINES record includes any of  the
 corresponding options.

 The following are equivalent:
 $SUBROUTINES ... TOL=n
 and
 $SUBROUTINES ...
 $TOL NRD=n

 Left-hand quantities:

 These are identical to the left-hand variables that may be set in sub-
 routine TOL.  If subscript "(I)" is omitted, it defaults to 1.

  NRD(I)
      The number of digits that are required to be accurate in the com-
      putation  of the drug amount in compartment I, i.e., the relative
      tolerance.  ADVAN 9, 13, 16, 17, and 18 have  the  capability  of
      using  different  values  of NRD for different compartments.  For
      compartments not specified, the tolerance of the last compartment
      specified will be used.

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
      10**(-12)).   Used  by  ADVAN  9,  13, 14, and 15, 16, 17, and 18
      which have the capability of using different values of  ANRD  for
      different compartments.  For compartments not specified, the tol-
      erance of the last compartment specified will be used.

      ANRD(0) is the absolute tolerance for the Steady  State  computa-
      tions.  If ANRD(0) is not specified, ANRD(1) is used.

      The  value  of  ANRD(1)  can also be specified using $SUBROUTINES
      option ATOL.
      The value of ANRD(0) can also  be  specified  using  $SUBROUTINES
      option SSATOL.

  NRDC(I) (NM74)
      Same  as  NRD(I),  but used for the FOCE/LAPLACE covariance step.
      Used with ADVAN 9, 13, 14, 15, 16, 17, and 18.  If not set,  NRDC
      defaults  to  the value of NRD.  NRDC(0) is used for Steady State
      computations during the FOCE/LAPLACE covariance step.

      The value of NRDC(1) can also  be  specified  using  $SUBROUTINES
      option TOLC.
      The  value  of  NRDC(0)  can also be specified using $SUBROUTINES
      option SSTOLC.

  ANRDC(I) (NM74)
      Same as ANRD(I), but used for the FOCE/LAPLACE  covariance  step.
      Used  with  ADVAN  9,  13, 14, 15, 16, 17, 18.  If not set, ANRDC
      defaults to the value of ANRD.  ANRDC(0) is used for Steady State
      computations during the FOCE/LAPLACE covariance step.

      The  value  of  ANRDC(1) can also be specified using $SUBROUTINES
      option ATOLC.
      The value of ANRDC(0) can also be  specified  using  $SUBROUTINES
      option SSATOLC.

 Right-hand quantities:

   Integers, e.g., NRD(1)=4

 Forbidden Variable Names:

   All other variables.

 RECORD ORDER:

 Follows $SUBROUTINES and $MODEL

 REFERENCES: Guide IV, section V.C.3 , V.C.10 
 REFERENCES: Guide VI, section VI.D , Figure 41
