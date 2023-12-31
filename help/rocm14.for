


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     PROBLEM ITERATION COUNTERS                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD_INT, ONLY: NPROB,IPROB,S1NUM=>IDXSUP(1),S2NUM=>IDXSUP(2), &
                           S1NIT=>NITR_SUP(1),S2NIT=>NITR_SUP(2), &
                           S1IT=>NITSUP(1),S2IT=>NITSUP(2)

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NPROB,IPROB,IDXSUP(2),NITR_SUP(2),NITSUP(2)

 DISCUSSION:

  NPROB
      The number of problems in the run.

  IPROB
      The number of the current problem.

      IPROB  is  1  unless  there are multiple problems (e.g., multiple
      $PROBLEM records in the NM-TRAN control stream).

  S1NUM
      The number of the active first-level superproblem.
      S1NUM is 0 if there is no active first-level superproblem.

  S2NUM
      The number of the active second-level superproblem.
      S2NUM is 0 if there is no active second-level superproblem.

  S1NIT
      The number of iterations for the active first-level superproblem.
      S1NIT is 0 if there is no active first-level superproblem.

  S2NIT
      The number of iterations for the active  second-level  superprob-
      lem.
      S2NIT is 0 if there is no active second-level superproblem.

  S1IT
      The  number  of  the  current iteration of the active first-level
      superproblem.
      S1IT is 0 if there is no active first-level superproblem.

  S2IT
      The number of the current iteration of  the  active  second-level
      superproblem.
      S2IT is 0 if there is no active second-level superproblem.

 These  variables  may  be used as right-hand quantities in abbreviated
 code for initialization/finalization blocks.

 (See $super).

 Location prior to NONMEM 7: rocm14

 REFERENCES: None.
