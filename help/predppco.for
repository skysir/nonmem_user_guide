


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           PREDPP MODULES                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Global variables in PREDPP
 CONTEXT: User-supplied routines

 DISCUSSION:

 These modules contain values that are (for the most part) communicated
 from PRED to its subroutines.  Prior to NONMEM  7  these  values  were
 COMMON  blocks.   Following is a list of the MODULE, the (old) COMMON,
 and a description of the variables.  (See NONMEM_modules).

 PREDPP read-only modules contain values  that  are  communicated  from
 PREDPP to various user-subroutines PK, ERROR, DES and AES.

 PROCM_INT (PROCM1)
      NEWIND, for PK and ERROR

 PROCM_REAL (PROCM2)
      Non-event dose time and derivatives, for PK

 PROCM_REAL (PROCM3)
      Initiating dose record (at non-event dose time), for PK

 PROCM_REAL (PROCM4)
      Compartment amounts and derivatives, for PK and ERROR

 PROCM_INT (PROCM5)
      Number and map of active etas, for PK and ERROR (and TRANS)

 PROCM_REAL (PROCM6)
      Theta vector from NONMEM, for DES and AES

 PROCM_REAL (PROCM7)
      EVTREC from NONMEM, for DES and AES

 PROCM_CHAR (PROCM8)
      Format statements, for all routines

 PROCM_REAL (PROCM9)
      Time at which the state-vector (in PROCM4) was last computed, for
      PK

 PROCM_INT (PROCMA)
      Indicator variables associated with model event  times,  for  all
      routines

 PROCM_INT (PROCMB)
      Flag  indicating  final call to DES or AES after advance to event
      or non-event time, for DES and AES

 PROCM_INT (PROCMC)
      Flag indicating compartment initialization call to PK, for PK

 The following are special modules containing values that are  communi-
 cated from user routines to PREDPP.

 PRCOM_INT (PRCOMG)
      Override  default settings in ADVAN6, ADVAN8, ADVAN9, ADVAN13,and
      SS6

 PRCOM_LOG (PRCOMU)
      Restore pre-1990 behavior of bio-availability fraction

 PRCOM_REAL (PRCOMW)
      Fudge factor for error tests in ADVAN7/SS7

 These modules contain values that are communicated from subroutines to
 PREDPP.

 PKERR_REAL (PRDPK1)
      Communicate  model  event  time  parameters computed by PK to the
      ERROR routine

 PRMOD_INT (PRDPK2)
      Flag indicating possible change in model event time parameters

 PRMOD_REAL (PRDPK3)
      Compartment initialization values from PK

 PRMOD_INT (PRDPK4)
      Values of ISSMOD and I_SS ("Initial Steady-State") flags.

 PROCM_INT (PRDDE1)
      Initialization values that are communicated from  subroutine  DES
      and AES to PREDPP.

 PRINFN (PRINFN)
      Contains  values of INFN-defined variables which are communicated
      to other PREDPP user-routines.

 Each has its own entry in the Help document.

 REFERENCES: Guide VI, section III.I , IV.D , Figure 14
