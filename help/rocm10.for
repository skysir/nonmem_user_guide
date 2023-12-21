


 +--------------------------------------------------------------------+
 |                                                                    |
 |                       SIMULATION: NREP,IREP                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables

 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_INT, ONLY: NREP,IREP=>NCREP

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: NREP,NCREP

 DISCUSSION:

 NREP The  number  of replications in the Simulation Step, given by the
      NSUBS option of the $SIMULATION record.

 IREP The number of the current replication.

 These variables may be used as right-hand  quantities  in  abbreviated
 code for initialization/finalization blocks.

 Location prior to NONMEM 7: rocm10

 REFERENCES: None.
