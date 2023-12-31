


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $PRIOR                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Provides instructions for the PRIOR subroutine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $PRIOR subroutine  [(conditional" "clause1 ), (conditional" "clause2 ) ... ]
        [DISPLAY[=ALL|CNT]]   [ICMAX=n]
        [argument1 , argument2 ...]

 SAMPLE:
 $PRIOR TNPRI (PROBLEM 2) PLEV=.9999 ISS=0 IVAR=1

 DISCUSSION:
 Optional. Specifies the use of the PRIOR feature of NONMEM.  Note that
 $PRIOR is a control record, not a block of  abbreviated  code.  There-
 fore,  only those options that are listed here may be used. E.g., ver-
 batim code may not be used. Options and arguments may be in any order,
 and may be on more than one line.

 OPTIONS:

 subroutine

   Required.  Either  TNPRI  or NWPRI.  The following options and argu-
   ments apply to calls to this subroutine.  Another subroutine  option
   (or another $PRIOR record) may follow, with a new set of options and
   arguments.  Use only NWPRI for the new $ESTIMATION methods of NONMEM
   7.

 conditional clause

   Optional.   One  or more conditions, within parentheses.  The condi-
   tions are "AND"ed together, i.e., the subroutine is called when  all
   conditions  in  the  conditional  clause are true.  If there is more
   than one conditional clause, the clauses are "OR"ed together,  i.e.,
   the  subroutine  is  called  if all the conditions in any one condi-
   tional clause are true. See the Examples, below.  Conditions may  be
   one or more of:

   ESTIMATION or SIMULATION
     Specifies  the  NONMEM  task  for  which  the  subroutine is to be
     called.  If omitted, PRIOR is called for all tasks (i.e., for  all
     values of ICALL).  ESTIMATION and SIMULATION may not both be spec-
     ified in the same conditional clause.  ESTIMATION may  be  spelled
     ESTIMATE or ESTM; SIMULATION may be spelled SIMULATE or SIML.

     May  also  be specified as ICALL=n, ICALL.EQ.n, or ICALLn, where n
     is 2 (ESTIMATION) or 4 (SIMULATION).

   PROBLEM=n
     Specifies the problem for which the subroutine is  to  be  called.
     May  also  be  specified as PROBLEM=n or PROBLEM.EQ.n or PROBLEMn.
     PROBLEM may also be coded as IPROB.  Instead of  =,  .EQ.  may  be
     used. Other permitted tests are .NE., .LT., .LE., .GT., and .GE.

 DISPLAY[=ALL|CNT]

   Optional.   The PRIOR subroutine will contain code to print items of
   interest in the NONMEM report.  This is to assist the user in check-
   ing that the $PRIOR record is working correctly.

   DISPLAY=ALL  is the default when only DISPLAY is present. Lines such
   as the following are printed with every call to PRIOR:

   PRIOR ICALL,IPROB,IREP,CNT:  2  2  0   -41.898951681785.

   With DISPLAY=CNT, lines such as the following are printed only  when
   PRIOR is called for simulation or estimation.

   PRIOR CNT:    -41.898951681785

 ICMAX=n

   Optional.   The  PRIOR  subroutine will set the given value in ICMAX
   prior to calling the subroutine.  (See PRIOR Simulation: ICMAX).

 Subroutine arguments

   Optional.  The arguments are described in the help entries for NWPRI
   and  TNPRI.  They must be coded excatly as shown, i.e., no abbrevia-
   tions.  Any argument that is omitted defaults to 0.

   ITYP, NSAM, ISS, PLEV, CNT

     Arguments for both NWPRI and TNPRI

   NTHETA, NETA, NEPS, NTHP, NETP, NEPP, NPEXP

     Arguments unique to NWPRI

   IFND, MODE, IVAR
     Arguments unique to for TNPRI

 EXAMPLES:

 $PRIOR TNPRI (ESTIMATION, PROB .GE.3 ) IFND=1

   TNPRI is called with the Estimation step of problems 3  and  higher.
   IFND is set to 1 with these calls; all other arguments are 0.

 $PRIOR TNPRI IFND=1
 (EST, PROB 3) (EST,PROB 4) (EST,PROB 5)

   Same  as  the preceding example, if the run consists of exactly five
   problems.  TNPRI is called with the Estimation step of problems 3, 4
   and  5.   IFND is set to 1 with these calls; all other arguments are
   0.

 $PRIOR TNPRI (EST, PROB 3) IFND=1
        TNPRI (EST, PROB 4) IFND=1
        TNPRI (EST, PROB 5) IFND=1

   Same as the preceding example.  TNPRI is called with the  Estimation
   step  of problem 3.  IFND=1 with this call.  Similarly with problems
   4 and 5.  Note that IFND must be specified independently  each  time
   to be 1.  With re-specification of the subroutine, all arguments are
   reset to 0.  This permits identical or different arguments with each
   usage of the subroutine.  NM-TRAN will warn if no argument is speci-
   fied in a subsequent specification of TNPRI or NWPRI,  in  case  the
   arguments  have been omitted by error.  One or more arguments may be
   set explicitly to 0 to prevent the warning.  When the subroutine  is
   specified  more  than  once,  then all specifications must be condi-
   tional, i.e., conditional clauses are requrired.

   Within the same run, one may use TNPRI with some tasks or  problems,
   and  NWPRI  with  other  tasks or problems.  There may be at most 10
   $PRIOR records per problem.

 (See TNPRI, TNPRI, prior).
 (See tnpri example,  nwpri example).

 REFERENCES: Guide Introduction_VI
