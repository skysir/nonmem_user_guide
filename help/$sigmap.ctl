


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          $SIGMAP,$SIGMAPD                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives prior information for sigmas
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SIGMAP  value1  [value2]  [value3] ...
 $SIGMAPD  value1  [value2]  [value3] ...

 SAMPLE:
 ; Prior to SIGMA (NEPPxNEPP=1x1 of them)
 $SIGMAP 0.05 FIX
 ; Set degrees of freedom of SIGMA Prior (one value per SIGMA block)
 $SIGMAPD (1 FIX)

 DISCUSSION:

 These  are  called the informative forms of the $SIGMA record, for use
 with the NWPRI utility.

 $SIGMAP gives prior information for elements of the SIGMA matrix.
 $SIGMAPD gives degrees of freedom (also called the dispersion  factor)
 for SIGMA priors.

 The  name  of  the  record describes the kind of information it gives,
 rather than the structure of the information.  E.g.,  in  the  example
 above, $SIGMAPD is implemented in FCON with a $THETA record because it
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
