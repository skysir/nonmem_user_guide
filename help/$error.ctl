


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $ERROR                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Marks the beginning of abbreviated code for the ERROR routine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $ERROR
 abbreviated code

 DISCUSSION:
 The  $ERROR record is used to model intra-individual error in observed
 values.  It is used with PREDPP.  It can also be used  to  to  convert
 predictions  from PREDPP, i.e., scaled drug amounts, to other types of
 predictions (for example, to obtain the prediction of a drug effect as
 a  function  of  concentration,  in a pharmacodynamic study).  General
 rules for abbreviated  code  are  documented  elsewhere  (See abbrevi-
 ated code).  Specific rules follow.

 ASSIGNMENT AND CONDITIONAL STATEMENTS

   Left-hand quantities in assignment statements:

     Y  (Required.  The  modeled value for the dependent variable under
     the statistical model.)

     ERROR-defined (i.e., PRED-defined) items.

     Left-hand quantities from the $PK block, if they  are  not  random
     variables and neither $PK nor $ERROR include the COMRES=-1 pseudo-
     assignment statement.

 Right-hand quantities in assignment statement and in conditions:

   F (Required. The value of the scaled drug amount in the  observation
   compartment.)

   Data item labels specified on the $INPUT statement.

   THETA(n)

   ETA(n)   (Required  if the data are single-subject, and can be coded
   ERR(n).  Optional if the data are population.)

   EPS(n)  (Required if the data  are  population,  and  can  be  coded
   ERR(n).)

   ERROR-defined  items  that appeared earlier as left-hand quantities.
   This includes Y.

   Left-hand quantities from the $PK block, if neither $PK  nor  $ERROR
   include the COMRES=-1 pseudo-assignment statement.

   NEWIND
    Same as the NEWIND argument passed by NONMEM to PREDPP.
    NEWIND=0:  First  record  of  the data set.  THETA value may differ
    from value at last call with this record.
    NEWIND=1: First record of the data set, THETA value does not differ
    from  value at last call with this record, and PRED is nonrecursive
    (see I_REC), or,
    First record of a subsequent individual record.
    NEWIND=2: Subsequent data record of an individual record.

   NEWL2
    NEWL2=1: First record of an L2 record.
    NEWL2=2: Otherwise.

   ICALL
    Same as the ICALL argument passed by NONMEM to PREDPP.
    ICALL=1: Initialization.
    ICALL=2: Normal call.
    ICALL=3: Finalization.
    ICALL=4: Simulation.
    ICALL=5: Expectation.
    ICALL=6: Data Average.
    Special rules apply to blocks of abbreviated code that are executed
    when ICALL is not 2.
    (See Initialization-Finalization block, Simulation block).
    (See Expectation block, Data_Average block).

 Global Variables in Modules
  Certain  variables  in  FORTRAN  Modules  can  be  used.   (See Vari-
  ables_in_modules) The following are of particular interest.

   A(n)  (Amount in compartment n.)  (See State Vector A)
    Note: If there is no verbatim code and no explicit use  of  a  sub-
    scripted  variable  A in the $ERROR block, then the symbol A can be
    used as a data item label or as a name of an ERROR-defined item.

 Forbidden Variable Names:
  IDEF IREV EVTREC NVNT INDXS G HH DADT(n) E(n) P(n)

 PSEUDO ASSIGNMENT STATEMENTS

  COMRES=-1

  CALLFL=-1: Call with every event record.
  CALLFL=0: Call with every observation record.
  CALLFL=1: Call once per individual record.

  Of the last three, CALLFL=-1 is the default, except when the abbrevi-
  ated  code  consists of one line specifying a simple additive or pro-
  portional error model, and no verbatim code is present.
  CALLFL=2: Call once per problem.
  In effect, this is supplied by NM-TRAN for the simple  error  models,
  but it cannot be specified explicitly in the $ERROR record.
  (See callfl).

  Pseudo-assignments statements may be enclosed in parentheses.  If two
  of them are present in the same set  of  parentheses,  separate  them
  with  a  semicolon.   A  calling  protocol  phrase may be used within
  parentheses instead of a pseudo-assignment statement, and then either
  upper or lower case may be used .  E.g.,
  $ERROR (ONCE PER IR)       ; same as CALLFL=1
  $ERROR (ONLY OBSERVATIONS) ; same as CALLFL=0
  $ERROR (EVERY EVENT)       ; same as CALLFL=-1 (default)

 RECORD ORDER:

   Follows $SUBROUTINES and $INPUT
   Follows $MODEL (if a general model such as ADVAN6 is used)
   Follows $PK (if present)

 REFERENCES: Guide IV, section V.C.6 
 REFERENCES: Guide V, section 8 
 REFERENCES: Guide VI, section IV 
