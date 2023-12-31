


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $ANNEAL                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING:  Sets  starting diagonal Omega values to facilitate EM search |
 methods
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $ANNEAL number-list1:value1 number-list2:value2 ...

 SAMPLE:
 $ANNEAL 1-3,5:0.3 6,7:1.0

 DISCUSSION:

 Sets starting diagonal Omega values for purposes of simulated  anneal-
 ing by NONMEM subroutine CONSTRAINT.

 In  the  above  example,  initial  values  of  OMEGA(1,1), OMEGA(2,2),
 OMEGA(3,3), and OMEGA(5,5) are set to 0.3,  while  initial  OMEGA(6,6)
 and OMEGA(7,7) are set to 1.0.

 A  number-list may contain a single integer, a range of integers (with
 -), or a series of integers and ranges separated by comma.  Required.

 A value may be any numeric value. Optional; default is 0.

 When $EST CONSTRAIN>=4, an algorithm  in  subroutine  CONSTRAINT  will
 initially  set the omegas to these values, and then shrink these OMEGA
 values more and more with each iteration, and eventually  shrinks  the
 OMEGA's  to  0,  the  intended target value for that Omega.  This is a
 technique that may be used especially with SAEM, to provide an anneal-
 ing  method for moving thetas that have 0 omega values associated with
 them.  The default is the use of gradient methods, which are good  for
 problems  starting  near the solution, whereas the annealing method is
 more suitable for problems starting far from the solution.

 Subroutine CONSTRAINT obtains values entered via $ANNEAL record in the
 array  OMEGANNL.   Any  value that is set to 0 in $ANNEAL is given the
 default value of .3.

 This record is optional. If omitted, the starting values of Omega  are
 those specified in the $OMEGA record.

 NONMEM's default CONSTRAINT.f90 is identical to Subroutine source/CON-
 STRAINT.f90.

 See INTRODUCTION TO NONMEM 7, $ANNEAL to facilitate EM search methods

 REFERENCES: Guide Introduction_7
