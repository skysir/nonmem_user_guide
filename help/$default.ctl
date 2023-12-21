


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $DEFAULT                              |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies certain defaults for NONMEM
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $DEFAULT  [NOSUB=[-1|0|1]]

 SAMPLE:
 $DEFAULT NOSUB=1

 DISCUSSION:
 $DEFAULT is optional.  If present, it must appear following $PROB

 Specifies  certain  defaults  for  NONMEM.   If more than one $DEFAULT
 record is present in a given problem, the one used by  NONMEM  is  the
 last one in the problem.

 OPTIONS:

 NOSUB=[-1|0|1]
      With  NOSUB=0, label substitution will be performed for all tasks
      in the problem.  This is the default.
      (See $ABBREVIATED).
      With NOSUB=1, label substitution will not be performed.
      With NOSUB=-1, revert to NONMEM default, which is to treat -1  as
      a 0.

      If  the  NOSUB  option  is also specified on a task specification
      record ($TABLE, $SCATTER), then this value of NOSUB applies  only
      for the current task.  When speficied on a $EST record, the usual
      rule for options apply, in which the option varies  carries  into
      the  next  $EST  record in the problem unless otherwise re-speci-
      fied.

 May also be coded $DEFAULTS.

 REFERENCES: Guide Introduction_7
