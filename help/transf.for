


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         TRANS (SUBROUTINE)                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: TRANS subroutine
 CONTEXT: User-supplied subroutine; for use with PREDPP

 USAGE:
 SUBROUTINE TRANS (ITRANS,IRGG,GG,NETAS)
 USE SIZES, ONLY: ISIZE,DPSIZE,LVR
 INTEGER(KIND=ISIZE) :: IRGG,NPETAS,ITRANS
 REAL(KIND=DPSIZE)   :: GG(IRGG,LVR+1,LVR+1)

 DISCUSSION:
 The  TRANS  subroutine translates (or transforms) the values for a set
 of basic PK parameters modeled in PK to a set of values for the inter-
 nal parameters required by the ADVAN.  The PREDPP Library has a number
 of TRANS subroutines, representing  different  possible  parameteriza-
 tions  in PK, from which the user may choose. If a suitable translator
 is not found in the Library, the user may write his own.

 Input/Output argument:

  ITRANS
      ITRANS=1: TRANS has been called for initialization at the  begin-
      ning of a NONMEM problem; one such call per problem.  ITRANS must
      be reset by TRANS to a number in the range 1-8999.   This  number
      appears on NONMEM output, allowing the user to identify the TRANS
      routine being used.

      ITRANS=2: On input the GG array has stored in it the values  com-
      puted  by  PK  (except  that were any PK parameter modeled in its
      logarithm form in PK, PREDPP would have already exponentiated its
      typical/subject-specific value and multiplied its eta derivatives
      by its exponentiated typical/subject-specific value).  On  output
      the  GG  array  should have stored in it the values that would be
      computed by PK were the internal parameters of the ADVAN  modeled
      directly in PK (and none in their logarithmic form).

      ITRANS=4: TRANS has been called during the Simulation Step.  Only
      the first column of  the  GG  array  need  be  computed  as  with
      ITRANS=2.  Other columns need not be computed.

  GG  The array of PK parameters and their eta derivatives.
      The array is described elsewhere; (See pk).

 Input argument:

  NETAS
      The number of population etas in the problem.

 Certain variables are available in modules.  Their use is optional.

 Variables in read/write modules:

   IERPRD NETEXT ETEXT
   (See PRED Exit Code).
   (See PRED Error Message).

 Variables in read-only commons:

   NOETAS, SECOND
   (See Partial Derivative Indicators).

 When  SECOND  is  true,  PREDPP requires second-partial derivatives of
 etas.

 REFERENCES: Guide VI, section III.M 
