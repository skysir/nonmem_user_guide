


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    PRIOR SIMULATION: PARAMETERS                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPR_REAL, ONLY: THSIMP=>THET_P,OMSIMP=>OMEG_P,SGSIMP=>SIGM_P,THSIMPR

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LTH,LVR,DPSIZE
      REAL(KIND=DPSIZE) :: THET_P(LTH),OMEG_P(LVR,LVR),SIGM_P(LVR,LVR)

 DISCUSSION:
 Values of THETA, OMEGA and SIGMA that are produced during a Simulation
 Step using the user-supplied routine PRIOR must  be  stored  in  these
 variables.   If  one  of the NONMEM utility routines NWPRI or TNPRI is
 used by PRIOR to produce THETA, OMEGA, and SIGMA values, these  values
 will  automatically  be  stored  in  THET_P,OMEG_P,SIGM_P.  If $THETAR |
 record is not  used,  THSIMPR  contains  the  same  values  as  THSIMP |
 (THET_P).   If $THETAR record, THSIMP contains the internal theta val- |
 ues, and THSIMPR contains the reported values  (THSIMP  values  trans- |
 formed  by  equations of the $THETAR record).  The entire THETA, OMEGA
 and SIGMA arrays are stored by NPWRI and TNPRI, not just the the  ini-
 tal sub-vectors (prior-affected parts).  (See nwpri, tnpri).

 These  variables  may be used as right-hand quantities in $PK, $ERROR,
 $INFN and $PRED blocks.  After being set during the  Simulation  Step,
 they remain available during problem finalization (i.e., ICALL=3).

 Location prior to NONMEM 7: nmpr16

 REFERENCES: None.
