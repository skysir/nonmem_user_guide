


 +--------------------------------------------------------------------+
 |                                                                    |
 |                 INITIAL STEADY STATE: I_SS,ISSMOD                  |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: PK routine

 USAGE:
      USE PRMOD_INT, ONLY: I_SS,ISSMOD

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: I_SS,ISSMOD

 DISCUSSION:

 Used  for the Initial Steady State feature of PREDPP, with the general
 non-linear models (ADVAN6, ADVAN8, ADVAN9, ADVAN13, ADVAN14,  ADVAN15,
 ADVAN16, ADVAN17, ADVAN18).

 By  default,  the  initial  conditions (i.e., compartment amounts) are
 zero at the start of each individual record.  Different initial condi-
 tions may be computed using the I_SS (Initial Steady State) feature of
 MODEL and/or PK.

 ISSMOD may be set by MODEL.
 (Note that the option of the $MODEL record is I_SS, not ISSMOD.)
 Default: -1 (MODEL does not request the I_SS feature)
 Values 0, 1, 2 and 3 are permitted.  Value 0 requests that  no  steady
 state be computed.  ISSMOD values 1, 2, or 3 requests that PREDPP com-
 pute an initial steady state for the  model  before  the  first  event
 record  of  an individual record, or after a reset event.  The results
 are identical to those that would be computed by a  steady-state  dose
 event  record with SS=ISSMOD and AMT=0 and RATE=0.  If endogenous drug
 is specified in the differential equations,  non-zero  initial  condi-
 tions will be computed.

 I_SS may be set by PK.
 Default: -1 (PK does not request the I_SS feature)
 Values are the same as ISSMOD, with the same effect.  This allows I_SS
 to be set conditionally, e.g., if some subjects  are  at  steady-state
 and  others  are not.  If both ISSMOD and I_SS are set, then the value
 of I_SS overrides that of ISSMOD.

 There is no difference between values 1, 2 and 3 of I_SS unless the PK
 routine  also  uses  the  compartment initialization feature A_0.  The
 I_SS feature behaves exactly like a steady state dose record  in  this
 regard.  Specifically,

      With I_SS=1 ("reset"), values of A_0 are ignored.
      With I_SS=2 ("sum"), values of A_0 are added to the SS values.
      With  I_SS=3  ("initial ests"), values of A_0 are used as initial
      estimates when computing the SS values.

 (and similarly for ISSMOD).

 If ISSMOD is set, or I_SS is set at ICALL=1,  PREDPP  will  acknowlege
 the fact in the NONMEM output.

 When I_SS or ISSMOD is set, even if only to 0, then SSS and the appro-
 priate SS routine are included in the NONMEM executable.

 When I_SS or ISSMOD is set, even if only to 0,  or  verbatim  code  is
 present  in  the $PK or $INFN block, then NM-TRAN includes I_SS,ISSMOD
 in the PK and/or INFN routine.

 Example: $PK may include code such as
   IF (TYPE.EQ.1) THEN
   I_SS=1
   ELSE
   I_SS=0
   ENDIF

 (See $model, $pk, model, pk).

 Location prior to NONMEM 7: prdpk4

 REFERENCES: None.
