


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                $AES                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Marks the beginning of abbreviated code for the AES routine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $AES
 abbreviated code

 DISCUSSION:
 The  $AES record is used to compute algebraic expressions which can be
 regarded as specifying equilibrium kinetics.  It is used with PREDPP's
 ADVAN9 and ADVAN15 and ADVAN17.
 (See $AESINITIAL).
 General rules for abbreviated code are documented elsewhere
 (See abbreviated code).
 Specific rules follow.

 ASSIGNMENT AND CONDITIONAL STATEMENTS

 Left-hand quantities in assignment statements:

   E(ncm1+1), E(ncm1+2), ...  (Required. Expressions which, when set to
   0, constitute the system of  algebraic  expressions  describing  the
   equilibrium  kinetics. ncm1 is the number of nonequilibrium compart-
   ments.)  The indices "(...)" may be omitted, in which  case  NM-TRAN
   will  supply  them  according  to the order in which the expressions
   appear.  Indices are required when the expressions are defined  con-
   ditionally (i.e., using an IF statement).

   AES-defined (i.e., PRED-defined) items.

 Right-hand quantities in assignment statement and in conditions:

   A(1),  A(2), ...   (Current compartment amounts; may be random vari-
   ables.)

   P(1), P(2), ...   (Post-translation basic PK parameters; may be ran-
   dom variables.)

   PK-defined  items (Implicit basic PK parameters; may be random vari-
   ables.)

   T (Time; may be random variable. T takes values continuously over an
   integration interval.)

   AES-defined  variables that appeared earlier as left-hand quantities
   in $AES, and  similarly  from  the  $AESINITIAL  record.   (Caution:
   AESINITIAL-defined variables that depend on compartment amounts will
   depend on the initial values of these compartment amounts,  not  the
   current values.)

   Data item labels specified on the $INPUT statement.

   THETA(n).

   Global Variables in Modules

   Certain variables in FORTRAN modules can be used.
   (See Variables_in_Modules)
   The following are of particular interest.

    DOSTIM
     DOSTIM  is  the  time of a lagged dose or additional dose to which
     the system is being advanced.  Abbreviated code in $AES  may  test
     DOSTIM.  It may use DOSTIM on the right, unless DOSTIM is a random
     variable.  However, it may be used on the right in a $PK block  to
     define a random variable which may in turn be used on the right in
     the $AES block.

    DOSREC
     DOSREC is the dose record corresponding to the  dose  entering  at
     DOSTIM.   Abbreviated  code  in $AES may test items in DOSREC in a
     logical condition, and DOSREC may always be used on the right.

    ISFINL
     During simulation or a copying pass, and during the advance  to  a
     particular  time  (event  or  non-event time), ISFINL=1 at a final
     call to AES at that time.  Otherwise, ISFINL=0.

 Forbidden Variable Names:

 IR DA DP DT ETA(n) EPS(n) ERR(n)

 RECORD ORDER:

 Follows $SUBROUTINES $INPUT $MODEL $PK
 Follows $AESINITIAL

 (See aes, advan9_15, advan9_17).

 REFERENCES: Guide IV, section V.C.9 
 REFERENCES: Guide VI, section VI.E 
