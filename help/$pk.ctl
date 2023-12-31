


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                $PK                                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Marks the beginning of abbreviated code for the PK Routine
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $PK
 abbreviated code

 DISCUSSION:
 The  $PK  record  is  used to model the values of basic and additional
 pharmacokinetic parameters.  It is used with PREDPP.

 Basic PK parameters are  typically  the  rate  constants  ("micro-con-
 stants") for use in kinetic formulas.  $PK can compute instead parame-
 ters such as clearance and volume, and a translator ("TRANS")  subrou-
 tine  can  be  used to convert these to rate constants.  Additional PK
 parameters include compartment scale parameters, which PREDPP uses  to
 convert compartment amounts to concentrations, and dose-related param-
 eters such as modeled infusion rates.  General rules  for  abbreviated
 code are documented elsewhere
 (See abbreviated code).
 Specific rules for $PK follow.

 ASSIGNMENT AND CONDITIONAL STATEMENTS

 Left-hand quantities in assignment statements:

   Basic  PK  parameters  for  the  ADVAN and TRANS routines (Required,
   except with the General Nonlinear Kinetics  Models  ADVAN6,  ADVAN8,
   ADVAN9,  ADVAN13, ADVAN14, ADVAN15, ADVAN16, ADVAN17, ADVAN18).  One
   or more of the following, depending on which  ADVAN  and  TRANS  are
   being used.
     K CL V (ADVAN1)
     KA K CL V (ADVAN2)
     K K12 K21 CL V Q VSS V1 V2 ALPHA BETA AOB (ADVAN3)
     KA K K23 K32 CL V Q VSS V1 V2 V3 ALPHA BETA AOB (ADVAN4)
     Km0 Kmn (ADVAN5, ADVAN7)
     P(n) (ADVAN6, ADVAN8, ADVAN9, ADVAN13, ADVAN14, ADVAN15, ADVAN16, ADVAN17, ADVAN18)
     VM, KM (ADVAN10)
     K K12 K21 K13 K31 CL Q2 Q3 V1 V2 V3 ALPHA BETA GAMMA
      (ADVAN11)
     KA K K23 K32 K24 K42 CL Q3 Q4 V2 V3 V4 ALPHA BETA GAMMA
      (ADVAN12)

     P(n) are referred to as "explicit" basic PK parameters.
     Any variable defined in $PK may be used on the right-hand side in a
     $DES or $AES block; these are "implicit" basic PK parameters.

   Additional PK parameters (Optional)
   One or more of the following,
   depending on the compartments defined for the model.
   The digit following the letter is the compartment number.
     Scale parameters Sn, e.g.: S1 S2 S3 S4 SC S0.
     Bio-availability fractions Fn, e.g.: F1 F2 F3.
     Output fractions Fn, e.g.: F2 F3 F4 F0 FO.
     Infusion rates Rn, e.g.: R1 R2 R3.
     Infusion durations Dn, e.g.: D1 D2 D3.
     Absorption lags ALAGn, e.g.: ALAG1 ALAG2 ALAG3.
     Time scale: TSCALE (may be written XSCALE).
     Model event times MTIME(i), e.g.: MTIME(1) MTIME(2).
     (The subscript i is not a compartment number.)
     (See MTIME)

   Initial compartment amounts (Optional), e.g..: A_0(1) A_0(2).
   (See Compartment Initialization Block)

   Initial steady state flag I_SS (Optional).
   (See advan68, advan9,  $model)
   (See Initial Steady State: I_SS,ISSMOD).

   PK-defined (i.e., PRED-defined) items

   Right-hand quantities in assignment statement and in conditions:

     Data item labels specified on the $INPUT statement.

     THETA(n).

     ETA(n) (Used if the data are population.)

     PK-defined items that appeared earlier as left-hand
     quantities.

     NEWIND
      Same as the NEWIND argument passed by NONMEM to PREDPP.
      NEWIND=0:
      First record of the data set.
      THETA value may differ from value at last call with this record.
      NEWIND=1:
      First record of the data set, THETA value does not differ from
      value at last call with this record, and PRED is nonrecursive
      (see I_REC), or,
      First record of a subsequent individual record.
      NEWIND=2:
      Subsequent data record of an individual record.

     ICALL
      Same as the ICALL argument passed by NONMEM to PREDPP.
      ICALL=1: Initialization.
      ICALL=2: Normal call.
      ICALL=3: Finalization.
      ICALL=4: Simulation.
      ICALL=5: Expectation.
      ICALL=6: Raw data averages.
      Special rules apply to blocks of abbreviated
      code that are executed when ICALL is not 2.
      (See Initialization-Finalization block, Simulation block).
      (See Expectation block, Data_Average block).

 Global Variables in Modules

  Certain variables in FORTRAN modules can be used.
  (See Variables_in_modules)
  The following are of particular interest.

  DOSTIM DOSREC(n)

   If PK is not being called at an additional or lagged dose time, then
   DOSTIM = 0 and all elements of DOSREC are 0.

   If PK is called at an additional or lagged dose time  t,  then  DOS-
   TIM=t
   (See Guide VI, Chapter III)
   DOSREC  contains a copy of the dose event record which initiated the
   additional  or  lagged  dose  (actually,  only  of  the  final  row:
   EVTREC(NVNT,*)).   Data  items  may be referred to by position or by
   label, e.g., DOSREC(1) or DOSREC(ID).  In DOSREC, TIME and all  user
   (concomitant)  data  items  have  values from the next event record.
   All other NONMEM/PREDPP reserved data items  have  values  from  the
   initiating  dose  event  record.   (The  $BIND record may be used to
   override this.)

  A(n) TSTATE

   A(n) are the latest computed compartment amounts, and TSTATE is  the
   time  at which they were computed.  That is, A(n) are the amounts at
   the previous event time, or if at a later time, but before the  time
   for which PK is being called, a lagged or additional dose was given,
   or a regular infusion was terminated, or a modeled  event  occurred,
   then  A(n)  are  the  amounts at the latest such time.  If there are
   population etas, and A(n) are used in the $PK abbreviated code, then
   any  $OMEGA  records  referring to etas explicitly used in this code
   should precede the $PK record, or if an $MSFI  record  is  used,  it
   should precede the $PK record and include the option NPOP=m.
   Note:  If  there  is  no verbatim code and no explicit use of a sub-
   scripted variable A in the $PK record, then the symbol A can be used
   as a data item label or as a name of a PK-defined item.

  A_0FLG

   A_0FLG  signals  a  record with which it is possible to initialize a
   compartment amount.
   (See Compartment Initialization Block)

 Forbidden Variable Names:

   IDEF IREV EVTREC NVNT INDXS IRGG GG NETAS DADT(n) E(n) EPS(n)

 PSEUDO ASSIGNMENT STATEMENTS

   COMRES=-1

   CALLFL=-2: Call with every event record and at additional and lagged
   dose times.
   CALLFL=-1: Call with every event record.
   CALLFL=0: Call with first event record and new TIME values.
   CALLFL=1: Call once per individual record.

   Of  the  last four, CALLFL=-2 is the default when DOSREC, DOSTIM, or
   MTIME are used  explicitly  in  the  abbreviated  code.   Otherwise,
   CALLFL=-1 is the default.

   The  pseudo-assignments  statements  may be enclosed in parentheses.
   If two of them are present in the same set of parentheses,  separate
   them with a semicolon.  A calling protocol phrase may be used within
   parentheses instead of a  pseudo-assignment  statement,  and  either
   upper or lower case may be used. E.g.,

   $PK (ONCE PER IR)      ; same as CALLFL=1
   $PK (NEW TIME)         ; same as CALLFL=0
   $PK (EVERY EVENT)      ; same as CALLFL=-1
   $PK (NONEVENT)         ; same as CALLFL=-2

 RECORD ORDER:
   Follows $SUBROUTINES and $INPUT
   Follows  $MODEL  (with  General  Nonlinear  Kinetics  Models ADVAN6,
   ADVAN8,  ADVAN9,  ADVAN13,  ADVAN14,  ADVAN15,   ADVAN16,   ADVAN17,
   ADVAN18).
   Precedes $ERROR (if present)

 REFERENCES: Guide IV, section V.C.5 
 REFERENCES: Guide V, section 7 
 REFERENCES: Guide VI, section III 
