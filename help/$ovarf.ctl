


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $OVARF                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies the weighting to the standard deviations of OMEGA
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $OVARF 0 (default)
 $OVARF [value ... ]
 $OVARF [(value)[xn]] ...
 SAMPLE:

 $OVARF 2.0 5.0

 Where  2.0 is specified for the first omega block, 5.0 for the second.
 If the corresponding $OLKJDF value is negative then this  is  argument
 STDSSP in user-defined OMEGA_STD_PRIORU.f90.

 DISCUSSION:

 The  $OVARF  is  a separate record that allows the user to specify the
 weighting (inverse variance) to the standard deviations LKJ decorrela-
 tion degrees of freedom for each OMEGA block.  Used with NUTS method.

 $EST OVARF over-rides $OVARF.

 REFERENCES: Guide Introduction_7
