


 +--------------------------------------------------------------------+
 |                                                                    |
 |                            $TTDF (NM75)                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies t-distribution degrees of freedom for theta
 CONTEXT: NM-TRAN Control Record

 USAGE:

 $TTDF 0 (default)
 $TTDF [value ... ]
 $TTDF [(value)[xn]] ...

 SAMPLE:
 $TTDF (2.0)x2 (3.0)x2
 sets  t-distribution degrees of freedom for thetas 1 and 2 to value of
 2, and degrees of freedom for thetas 3 and 4 to value of 3.0.

 DISCUSSION:
 The ttdf is a separate record that allows the user to specify  the  t-
 distribution  degrees  of freedom for each theta.  It can be used with
 $ESTIMATION METHOD=NUTS, and with $SIMULATION.  By default $TTDF value
 is  0.   These values will be used in estimation and simulation unless
 $EST TTDF or $SIM TTDF is set, respectively, which set the  t  degrees
 of freedom to a single TTDF value for all thetas for that step.

 REFERENCES: Guide IV, section V.C.3 , V.C.10 
 REFERENCES: Guide VI, section VI.D , Figure 41
