


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          $OMEGAP,$OMEGAPD                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives prior information for omegas
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $OMEGAP  value1  [value2]  [value3] ...
 $OMEGAPD  value1  [value2]  [value3] ...

 SAMPLE:
 ; Prior to OMEGA (NETPxNETP=4x4 of them)
 $OMEGAP BLOCK(4)
 0.2 FIX
 0.0  0.2
 0.0  0.0 0.2
 0.0  0.0 0.0 0.2
 ; Set degrees of freedom of OMEGA Prior (one value per OMEGA block)
 $OMEGAPD (4 FIX)

 DISCUSSION:

 These  are  called  the  informative  forms  (i.e., informative record
 names) of the $OMEGA record, for use with the NWPRI utility.

 $OMEGAP gives prior information for elements of the OMEGA matrix.
 $OMEGAPD gives degrees of freedom (also called the dispersion  factor)
 for OMEGA priors.

 The  name  of  the  record describes the kind of information it gives,
 rather than the structure of the information.  E.g.,  in  the  example
 above, $OMEGAPD is implemented in FCON with a $THETA record because it
 is always a vector of values.  These records may be  located  anywhere
 in  the  control stream.  NM-TRAN inserts the corresponding records in
 the control stream in the correct order.  When the  informative  forms
 are used, the options of $PRIOR NWPRI need not be specified.  However,
 if options are listed in $PRIOR NWPRI, then these  values  are  chosen
 over what is surmised from the informatively labeled theta/omega/sigma
 records.
 (See nwpri).

 OPTIONS:

 The option FIXED should be used.  Other appropriate options are  BLOCK
 and VALUES.

 REFERENCES: Guide Introduction_7
