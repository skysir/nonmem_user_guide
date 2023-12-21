


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        EVENT RECORD: EVTREC                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied DES and AES routines

 USAGE:
      USE PROCM_REAL, ONLY: EVTREC
      USE PROCM_INT, ONLY: NVNT=>NEVENT

 GLOBAL DECLARATION:
      USE SIZES, ONLY: PD,DPSIZE
      REAL(KIND=DPSIZE) :: EVTREC(5,PD+1)
      INTEGER(KIND=ISIZE) :: NEVENT

 DISCUSSION:

 EVTREC
      Identical to the argument EVTREC passed by PREDPP to PK and ERROR
      routine

 NVNT Identical to the argument NVNT passed by PREDPP to PK  and  ERROR
      routines.

 Location prior to NONMEM 7: procm7

 REFERENCES: Guide VI, section III.C 
