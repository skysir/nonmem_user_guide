


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                SKIP                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPR_INT, ONLY SKIP_=>SKIP

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: SKIP

 DISCUSSION:

 The  SKIP_  variable controls premature termination of a problem (with
 subproblems), superproblem or superproblem iteration.  It may  be  set
 as follows during (sub)problem-finalization.

 SKIP_=1
      With  a  subproblem:  terminate the entire problem and proceed to
      the next problem, if this exists.

 SKIP_=2
      With a problem during an iteration of a second-level superproblem
      A nested within a first-level superproblem B: terminate the iter-
      ation and proceed to the next iteration, if this exists,  and  if
      this does not exist, to the following (second-level super)problem
      of B, if this exists, and if this does not  exist,  to  the  next
      iteration of B.

 SKIP_=3
      With a problem during an iteration of a second-level superproblem
      A nested within  a  first-level  superproblem  B:  terminate  the
      entire  superproblem A and proceed to the following (second-level
      super)problem of B, if this exists, and if this does  not  exist,
      to the next iteration of B.

 SKIP_=4
      With a problem during an iteration of a first-level superproblem:
      terminate the iteration and proceed to  the  next  iteration,  if
      this exists, and if this does not exist, to the following (first-
      level super)problem, if this exists.

 SKIP_=5
      With a problem during an iteration of a first-level superproblem:
      terminate  the  entire  superproblem and proceed to the following
      (first-level super)problem, if this exists.

 With NM-TRAN, in a finalization block of abbreviated code one may  set
 SKIP_ and/or use the following phrases:

 END PROBLEM                              (same as SKIP_=1)
 END SECOND-LEVEL SUPERPROBLEM            (same as SKIP_=3)
 END FIRST-LEVEL SUPERPROBLEM             (same as SKIP_=5)
 END SECOND-LEVEL SUPERPROBLEM ITERATION  (same as SKIP_=2)
 END FIRST-LEVEL SUPERPROBLEM ITERATION   (same as SKIP_=4)

 E.g. The following are all equivalent:

 A.
      SKIP_=2

 B.
      IF (...) THEN
         END SECOND-LEVEL SUPERPROBLEM ITERATION
      ENDIF

 C.
      IF (...) THEN
         ENDSECONDLEVELSUPERPROBLEMITERATION
      ENDIF

 Location prior to NONMEM 7: nmpr15

 REFERENCES: None.
