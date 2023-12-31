


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $MODEL                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies the MODEL subroutine of PREDPP
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $MODEL   [NCOMPARTMENTS=n1] [NEQUILIBRIUM=n2] [NPARAMETERS=n3]
          [COMPARTMENT=([name] [attribute1] [attribute2] ...)] ...
          [LINK compnamea [TO|AND] compnameb BY k [l]] ...
          [I_SS=n]

 SAMPLE:
 $MODEL   NPARAMETERS=3
          COMP=(DEPOT DEFDOSE INITIALOFF) COMP=(CENTRAL DEFOBS NOOFF)

 DISCUSSION:
 Required with a general ADVAN (ADVAN5,6,7,8,9,13,14,15,16,17,18).

 OPTIONS:

 NCOMPARTMENTS=n1
      Total  number  of compartments other than the output compartment.
      Default: the number of COMPARTMENT options.  May  also  be  coded
      NCM or NCOMPS.

 NEQUILIBRIUM=n2
      Number of equilibrium compartments.  Default: 0

 NPARAMETERS=n3
      Number  of  basic PK parameters.  Default: The number of basic PK
      parameters defined in the $PK  abbreviated  code.   May  also  be
      coded  NPARAMS.   May  be  0  with  the general non-linear models
      (ADVAN6, ADVAN8,  ADVAN9,  ADVAN13,  ADVAN14,  ADVAN15,  ADVAN16,
      ADVAN17, ADVAN18).

 COMPARTMENT= ([name] [attribute1] [attribute2] ...)
      Each  COMPARTMENT  option defines a single compartment.  Compart-
      ments are numbered in the order in which they are  defined.   The
      name  option gives the name for the compartment.  The name option
      may not be one of the compartment attributes below, unless it  is
      enclosed in quotes (' or ").
      If  the  name  option  is not used, then the compartment is named
      "COMPn," where n is the compartment number.

      Attributes are chosen from the following list.  When an attribute
      is not chosen, its opposite is the default.

        INITIALOFF - Compartment is initially off.
        NOOFF - Compartment may not be turned on or off.
        NODOSE - Compartment may not receive a dose.
        EQUILIBRIUM   -   Compartment  is  an  equilibrium  compartment
        (ADVAN9, ADVAN15, and ADVAN17 only; implies NODOSE).
        EXCLUDE - Compartment amount is excluded from  the  computation
        of the amount (Amt) of the output compartment (ADVAN9, ADVAN15,
        and ADVAN17 only).  [Amt=(total amount of drug thus  far  input
        into system) minus (total amount remaining in the system)]

      The following two attributes have no opposites.

        DEFOBSERVATION  -  Compartment  is the default observation com-
        partment.
        DEFDOSE - Compartment is the default dose compartment.
        If no compartment has the DEFOBSERVATION attribute, the default
        is  the first compartment with the name CENTRAL; otherwise, the
        first  compartment.   If  no  compartment   has   the   DEFDOSE
        attribute,  the  default is the first compartment with the name
        DEPOT which may receive doses; otherwise, the first compartment
        which may receive doses.

 LINK compname [TO|AND] compname BY k [l]
      Link  clauses  are  used  with ADVAN5 and ADVAN7, but only in the
      absence of $PK abbreviated code.  They are rarely used,  and  are
      described in Guide IV.

 I_SS= n
      I_SS  may  be  used  with  the general non-linear models (ADVAN6,
      ADVAN8, ADVAN9, ADVAN13, ADVAN14, ADVAN15 ) to request  the  ini-
      tial state feature of PREDPP.  Values of n are
        0 No initial state state (the default)
        1 Initial steady state
        2 Initial steady state, adds to current compartment amounts.
        3 Initial steady state, use current compartment amounts as ini-
        tial estimates.

 The attributes INITIALOFF NODOSE define  an  output-type  compartment.
 As with the output compartment, it is possible to use a negative value
 of CMT to obtain an observation and turn off such a  compartment  with
 one  observation  record. Like the output compartment, it is initially
 off, and it remains off (so that  the  amount  therein  remains  zero)
 until  it  is explicitly turned on by an other type event record which
 has the  output compartment's number in the CMT data item. Unlike  the
 output  compartment,  a  differential equation (DADT) must be supplied
 for an output-type compartment.

 (See advan68, advan9_15, $pk)
 (See Initial Steady State: I_SS,ISSMOD).

 RECORD ORDER:
 Follows $SUBROUTINES
 Precedes $PK (if present)

 SYMBOLIC LABEL SUBSTITUTIONS OF MODEL COMPARTMENTS  (NM75)

 Compartment names defined in $MODEL are  automatically  available  for
 substitution  in  abbreviated code.  This is referred to as "implicit"
 compartment name substitution.  It is  an  alternative  to  "explicit"
 compartmet name substitution using $ABBR REPLACE records

 For example:

 $MODEL
 COMP=(GUT,DEFDOS)
 COMP=(CENTRAL,DEFOBS)
 COMP=(PERI)

 allows  susbsitutions  to  be  made for A(GUT), DADT(GUT), A(CENTRAL),
 DADT(CENTRAL), etc, so you may use these symbols in  your  abbreviated
 code, as in the following:

 $DES
 DADT(GUT)=-KA*A(GUT)
 DADT(CENTRAL)=KA*A(GUT)-(KCP+KC0)*A(CENTRAL)+KPC*A(PERI)
 DADT(PERI)=KCP*A(CENTRAL)-KPC*A(PERI)
  ...
 $ERROR
 IPRED=A(CENTRAL)/S2

 In the generated code, compartment names are replaced by the appropri-
 ate values 1, 2, etc.

 See INTRODUCTION TO NONMEM 7, Symbolic Label  Substitutions  of  Model
 Compartments

 REFERENCES: Guide IV, section V.C.4 
 REFERENCES: Guide VI, section VI.B 
