


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                $MIX                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Marks the beginning of abbreviated code for the MIX routine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $MIX
 abbreviated code

 SAMPLE:

 $MIX
 NSPOP=2
 P(1)=THETA(5)
 P(2)=1-THETA(5)

 DISCUSSION:
 The  $MIX  record  is  used to describe the mixture probabilities of a
 mixture model. It is evaluated with individual records.

 General  rules  for  abbreviated   code   are   documented   elsewhere
 (See abbreviated code). Specific rules follow.

 ASSIGNMENT AND CONDITIONAL STATEMENTS

 Left-hand quantities in assignment statements:

   NSPOP

   An integer variable.  The number of sub-populations. Must be given a
   value when ICALL=1.

   P(i)

   For each i (i=1, ... , NSPOP), P(i) is the modeled fraction  of  the
   population  in  the  ith  subpopulation.  The sum of the P(i) should
   equal 1. In principle, the P(i) can change from individual to  indi-
   vidual.  If for a given individual, the second (for example) subpop-
   ulation doesn't apply, then set P(2)=0 for  that  individual.   MIXP
   may be coded instead of P; they are synonymous in the context of the
   $MIX block.

 Right-hand quantities in assignment statement and in conditions:

   THETA(n).

   MIX-defined items that appeared earlier as left-hand quantities.

   ICALL
     ICALL=0: Run initialization.
     ICALL=1: Problem initialization.
     ICALL=2: Normal call.
     ICALL=4: Simulation.
   Data items listed in DATA option of $CONTR record
     E.g., assume that the following  records  appear  in  the  control
     stream prior to the $MIX block:

     $INPUT ... STDY ...
     $CONTR DATA=(STDY)
     Then STDY may be used on the right in $MIX.  STDY and STDY(1) both
     refer to the value of STDY on the first observation record of  the
     individual  record.  STDY(i) refers to the value of STDY on the i-
     th. observation record of the individual record.   TEMPLT  may  be
     used.
     (See MIX CONTR: TEMPLT).

 Global Variables in Modules
  Certain  variables  in  FORTRAN  modules  can  be  used.   (See Vari-
  ables_in_modules)

 May not include:

   EXIT, CALL , DO WHILE / ENDDO statements

   COMRES, CALLFL pseudo-variables

 Forbidden Variable Names:

   NEWIND,ETA(i),EPS(i),ERR(i) COM(i)

 Variables defined in $MIX are not listed in module NMPRD4 and may  not
 be displayed in $TABLE and $SCATTER.

 RECORD ORDER:

 Follows $SUBROUTINES.

 REFERENCES: Guide IV, section III.B.4 , III.B.6 
 REFERENCES: Guide IV, section IV.E.1 , 4.E.2 
 REFERENCES: Guide VI, section III.L.2 , Figure 6
