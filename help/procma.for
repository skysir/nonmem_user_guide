


 +--------------------------------------------------------------------+
 |                                                                    |
 |                 MODEL EVENT TIME: MNOW,MPAST,MNEXT                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_INT, ONLY: MNOW=>MTNOW,MPAST=>MTPAST,MNEXT=>MTNEXT

 GLOBAL DECLARATION:
      USE SIZES, ONLY: PCT
      INTEGER(KIND=ISIZE) :: MTNOW,MTPAST(PCT),MTNEXT(PCT)

 DISCUSSION:

 These  variables  are of interest when model event time parameters are
 defined in PK.

 (See mtime).

 Location prior to NONMEM 7: procma

 REFERENCES: None.
