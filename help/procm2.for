


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    DOSE TIME NON-EVENT: DOSTIM                     |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied PK, DES, AES routines

 USAGE:
      USE PROCM_REAL, ONLY: DOSTIM,DDOST,D2DOST

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR,DPSIZE
      REAL(KIND=DPSIZE) :: DOSTIM,DDOST(LVR),D2DOST(LVR,LVR)

 DISCUSSION:

 DOSTIM
      DOSTIM=0:  this  call  to  PK  occurs  at  an event time.  DDOST,
      D2DOST, and DOSREC contain zeros.

      DOSTIM>0: this call to PK occurs at a non-event dose time DOSTIM,
      i.e.,  at  the  time  of  an  additional  or lagged dose.  DDOST,
      D2DOST, and DOSREC contain values of interest.  (See Dose Record:
      DOSREC).

 DDOST
      DDOST(i) = Partial derivative of DOSTIM with respect to eta(i).

 D2DOST
      D2DOST(i,j) = Second partial derivative of DOSTIM with respect to
      eta(i), eta(j) (lower-triangular; j=1, ..., i).

 When DOSTIM>0, TIME and all user (concomitant) data  items  in  EVTREC
 are from the next event record.  All other NONMEM/PREDPP reserved data
 items are from the initiating dose event record.   (The  $BIND  record
 may be used to override this; (See bind)

 Location prior to NONMEM 7: procm2

 REFERENCES: Guide VI, section III.I 
