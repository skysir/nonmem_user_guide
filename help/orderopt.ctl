


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            ORDER_OPTION                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Instructions for the NONMEM Estimation Step
 CONTEXT: NM-TRAN Control Record

 ORDER is an option of the $ESTIMATION record:

 ORDER=xxxf
      The  values of x may be T (Theta), S (Sigma), and O (Order).  The
      value of f may be U (Upper) or L (Lower).  Affects the way theta,
      omega,  and  sigma are displayed in the raw and additional output
      files.  xxx gives the overall order, and f gives the order within
      OMEGA  and SIGMA.  Affects the raw output file and all additional
      output files.  The default is TSOL: THETA, SIGMA, OMEGA in  Lower
      triangular form.  Does not affect the NONMEM report file.

 DISCUSION:

 On  the $OMEGA record, elements of Omega are given in lower triangular
 order, e.g.,
 $OMEGA BLOCK(3) OM11 OM21 OM22 OM31 O32 OM33

 This can also be coded as
 $OMEGA BLOCK(3)
 OM11
 OM21 OM22
 OM31 O32 OM33

 The NONMEM report file is  not  consistent.   The  initial  and  final
 parameter estimates are printed in Lower triangular order, e.g.,
 INITIAL ESTIMATE OF OMEGA:
 OM11
 OM21 OM22
 OM31 OM32 OM33

 However,  in  the  intermediate printout of the Estimation step, OMEGA
 and SIGMA are in Upper triangular  order.  This was not  evident  with
 previous  version of NONMEM, because the PARAMTER values corresponding
 to OMEGA and SIGMA are displayed as  unconstrained  parameters  (UCP),
 and  these  are  not  in  a 1-1 mapping with the elements of OMEGA and
 SIGMA.  But with NONMEM 7.2 and higher, the intermediate output of the
 Estimation  step  includes  a line NPARAMETR, in which OMEGA and SIGMA
 are converted to their "natural" space.  E.g.,
 OM11 OM12 OM13 OM22 OM23 OM33
 Because OMEGA is symmetric, this is the same as
 OM11 OM21 OM31 OM22 OM32 OM33

 In the NONMEM report, upper triangular form is used  for  all  of  the
 "COVARIANCE  MATRIX OF THE ESTIMATE", "CORRELATION MATRIX OF ESTIMATE"
 and the "INVERSE COVARIANCE MATRIX OF ESTIMATE".  E.g.,
 OM11 OM12 OM13 OM22 OM23 OM33

 This is TOSU order.

 Case 1. ORDER is not used (Default)
      The raw output file .ext and additional output files  .cov,  etc,
      are in a different order: TSOL.  SIGMA precedes OMEGA.  The order
      within  OMEGA is lower triangular, consistent  with  the  Initial
      and final estimates of OMEGA, but different than the intermediate
      output in the NONMEM report.  The order of the last line  of  the
      "COVARIANCE MATRIX OF ESTIMATE" in the NONMEM report differs from
      the order of the lines of the .cov file, and similarly  for  .cor
      and .coi.

 Case 2. ORDER=TOSU
      Each  line of the raw output file .ext is in the same order as in
      the intermediate output portion of the NONMEM report.  Also, each
      line  of  the  .cov file is in the same order as the last line of
      the "COVARIANCE MATRIX OF ESTIMATE" in  the  NONMEM  report,  and
      similarly for .cor and .coi.

 Note  that  the  ORDER  option does not affect the NONMEM report file,
 only the raw and additional output files.   The  SIGMA  matrix  always
 printed in the same order as OMEGA, and the ORDER option applies to it
 as well.

 REFERENCES: Guide Introduction_7
