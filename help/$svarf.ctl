


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $SVARF                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies the weighting to the standard deviations of SIGMA
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SVARF 0 (default)
 $SVARF [value ... ]
 $SVARF [(value)[xn]] ...
 SAMPLE:

 $SVARF 2.0 5.0

 Where  2.0 is specified for the first sigma block, 5.0 for the second.
 If the corresponding $SLKJDF value is negative then this  is  argument
 STDSSP in user-defined SIGMA_STD_PRIORU.f90.

 DISCUSSION:

 The  $SVARF  is  a separate record that allows the user to specify the
 weighting (inverse variance) to the standard deviations LKJ decorrela-
 tion degrees of freedom for each SIGMA block.  Used with NUTS method.

 $EST SVARF over-rides $SVARF.

 REFERENCES: Guide Introduction_7
