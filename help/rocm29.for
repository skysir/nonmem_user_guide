


 +--------------------------------------------------------------------+
 |                                                                    |
 |                POPULATION SINGLE-SUBJECT INDICATOR                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD_INT, ONLY: IPS

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IPS

 DISCUSSION:

  IPS IPS = 1 when the data are population data.
      IPS = 2 when the data are single-subject data.

 Location prior to NONMEM 7: rocm29

 REFERENCES: Guide I, section B.1 
 REFERENCES: Guide IV, section II.C.4 
