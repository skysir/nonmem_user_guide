


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          $THETAP,$THETAPV                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives prior information for thetas
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $THETAP  value1  [value2]  [value3] ...
 $THETAPV  value1  [value2]  [value3] ...

 SAMPLE:
 ; Prior information of THETAS (NTHP=4 of them)
 $THETAP (2.0 FIX) (2.0 FIX) (2.0 FIX) (2.0 FIX)
 ; Variance to prior information of THETAS (NTHPxNTHP=4x4 of them).
 $THETAPV BLOCK(4)
 10000 FIX
 0.00 10000
 0.00  0.00 10000
 0.00  0.00 0.0 10000

 DISCUSSION:

 These  are  called the informative forms of the $THETA record, for use
 with the NWPRI utility.

 $THETAP gives prior information for elements of the THETA matrix.
 $THETAPV gives variance information for THETA priors.

 The name of the record describes the kind  of  information  it  gives,
 rather  than  the  structure of the information.  E.g., in the example
 above, $THETAPV is implemented in FCON with an $OMEGA  record  because
 it  is  in  general  an array of values.  These records may be located
 anywhere in the control stream.   NM-TRAN  inserts  the  corresponding
 records in the control stream in the correct order.  When the informa-
 tive forms are used, the options of $PRIOR NWPRI need  not  be  speci-
 fied.  However, if options are listed in $PRIOR NWPRI, then these val-
 ues are chosen over what is surmised from  the  informatively  labeled
 theta/omega/sigma records.
 (See nwpri).

 OPTIONS:

 The  option FIXED should be used.  Other appropriate options are BLOCK
 and VALUES.

 REFERENCES: Guide Introduction_7
