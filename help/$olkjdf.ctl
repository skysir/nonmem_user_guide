


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $OLKJDF                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies LKJ decorrelation degrees of freedom for each OMEGA
 block
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $OLKJDF 0 (default)
 $OLKJDF [value ... ]
 $OLKJDF [(value)[xn]] ...

 DISCUSSION:
 OLKJDF is an option of the $ESTIMATION record.  $OLKJDF is a  separate
 record  that  allows  the user to specify LKJ decorrelation degrees of
 freedom for each OMEGA block.

 SAMPLE:
 $OLKJDF 4.5 3.5 -2.0
 Where 4.5 degrees of freedom are specified for the first omega  block,
 3.5  for  the second, and - 2.0 specifies 2 degrees of freedom for the
 third omega block, but a user-defined definition of the standard devi-
 ations   of   the  diagonals  for  the  third  omega  block.  Use  the
 OMEGA_STD_PRIORU.f90 file in the ..\source directory as a template  to
 modify the OMEGA_STD_PRIORU to provide the desired probability density
 function. Set $SUBR OTHER=OMEGA_STD_PRIORU.f90 to use the  user  modi-
 fied template (you may rename the file for organizational purposes).

 RECORD ORDER:

 Follows $OMEGA

 REFERENCES: Guide Introduction_7
