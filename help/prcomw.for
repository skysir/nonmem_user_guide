


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           ADVAN7: FACTOR                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP global variables
 CONTEXT: For use with PREDPP

 USAGE:
      $PK
      "FIRST
      " USE PRCOM_REAL, ONLY: FACTOR=>FAC
      "MAIN
      " FACTOR=10.

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: FAC

 DISCUSSION:

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 PREDPP code.

 USAGE:

 A message similar to the following is sometimes produced by ADVAN7:

 NUMERICAL DIFFICULTIES OBTAINING THE SOLUTION.  NON-REAL
 EIGENVALUES. PERHAPS ADVAN5 SHOULD BE USED FOR THIS PROBLEM.
 "LARGEST" VALUE ALONG SUBDIAGONAL OF SCHUR MATRIX:  0.423E-01
 (DIAG: 0.84E+00 0.84E+00)

 FACTOR  is a "fudge factor" that affects ADVAN7's determination of how
 large the subdiagonal element is relative to  the  two  diagonal  ele-
 ments.  The default value of FACTOR is 1.

 These  error  messages  arise when the computations of the eigenvalues
 produce unexpected values.  The question is  whether  the  subdiagonal
 elements are too large to ignore (e.g., when the eigenvalues are truly
 complex) in which case ADVAN5 must be used, or whether trivially small
 quantities have been produced by the accumulation of rounding error in
 the floating point computation, in which ADVAN7 can  be  used.   (Non-
 real eigenvalues result when the system is non-mamillary, and may also
 occasionally occur when there are etas associated with rate  constants
 Kij).

 If the maximum value printed is very small (e.g., E-15), then the com-
 putations will proceed correctly using a value of FACTOR such as 10 or
 100.   Values such as FACTOR=10 or FACTOR=100 are more permissive, and
 allow larger subdiagonal elements.  If no reasonable value  of  FACTOR
 eliminates the error messages, then ADVAN5 should be used.

 Location prior to NONMEM 7: prcomw

 REFERENCES: None.
