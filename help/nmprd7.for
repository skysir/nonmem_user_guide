


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        SIMULATION: ETA,EPS                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_REAL, ONLY: ETA,EPS

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR,DPSIZE
      REAL(KIND=DPSIZE) :: ETA(LVR),EPS(LVR)

 DISCUSSION:

 When  the ONLYSIMULATION option on the $SIMULATION record is used, the
 displayed eta values are those that PRED stores in ETA during the Sim-
 ulation Step (i.e., at ICALL=4).

 In  generated subroutines, eta values are automatically stored in ETA,
 and thus automatically, the simulated etas are displayable.

 When using the PRED repetition feature, epsilon  values  generated  in
 PRED for simulation purposes (with or without the use of the ONLYSIMU-
 LATION option) should be  stored  in  EPS.    Then  when  records  are
 repeated,  these  same  epsilons  values  will  be available in EPS as
 input.  In particular, with every repeated  record,  the  values  that
 were  stored  in  EPS  the  last  time the record was passed, are made
 available as input to PRED.

 In generated subroutines, epsilon values are automatically  stored  in
 EPS,  and  thus automatically, the simulated epsilons are available as
 input with repeated records.

 When using the repetition feature with single-subject data, eta values
 generated  in  PRED  for  simulation purposes should be stored in ETA.
 What happens in this case is analogous to what happens in the case  of
 population data with epsilon values.

 In  generated subroutines, eta values are automatically stored in ETA,
 and thus automatically, the simulated etas are available as input with
 repeated records.

 Location prior to NONMEM 7: nmprd7

 REFERENCES: None.
