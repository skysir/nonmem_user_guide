


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    NON-ACTIVE ETA LIST FOR PRED                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY NAETA,LVOUT

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LVR
      INTEGER(KIND=ISIZE) :: NAETA,LVOUT(LVR)

 DISCUSSION:

 NAETA
      The number of positions in LVOUT that are set to 1.

 LVOUT
      When  LVOUT(i)=1,  NONMEM  is  ignoring  partial derivatives with
      respect to eta(i) with the current call to PRED.

 Location prior to NONMEM 7: rocm15

 REFERENCES: None.
