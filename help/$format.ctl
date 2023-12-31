


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $FORMAT                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies significant digits for the NONMEM report file
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $FORMAT  FMTN=n (NM75)

 $FORMAT FMTN=3 (default) (NM75)

      You  may  now  display any significant digits for thetas, omegas,
      and sigmas results, and variance- covariance of estimates, listed
      in  the NONMEM report file. The FMTN is the number of significant
      digits (which is by legacy and default 3), between 3 and 23. This
      also  applies  to  $TABLE  outputs to the NONMEM report file. The
      FORTRAN format will be formed as 1PE{FMTN+6}.{FMTN-1}. For  exam-
      ple  FMTN=6  will be 1PE12.5. If you wish G field format, set the
      FMTN to a negative value. For example FMTN=-4 will be 1PG10.3.

      It is recommended that you place the $FORMAT  record  immediately
      after the $PROB record:

      $PROB My problem
      $FORMAT FMTN=5

      This  FMTN  precision format will be in effect for outputs of all
      problems, until another $FORMAT record is given for a  new  $PROB
      in the control stream.

 REFERENCES: Guide Introduction_7
