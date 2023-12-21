


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        MIXTURE MODEL: MIXP                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_REAL, ONLY: MIXP

 GLOBAL DECLARATION:
      USE SIZES, ONLY: MMX,DPSIZE
      REAL(KIND=DPSIZE) :: MIXP(MMX)

 DISCUSSION:

 Used with mixture models.

  MIXP(i)
      The  mixture  probability  for the ith subpopulation, computed by
      subroutine MIX.

 May be used as right-hand quantities in all  abbreviated  codes  other
 than $MIX.  There the index i must be an integer constant.

 This is a general description. For details, see mixnum.                |

 (See MIXNUM_MIXEST_MIXP)
 (See mixnum_mixest_for_mixture_model)
 (See mix, mixture model example)

 Location prior to NONMEM 7: rocm25

 REFERENCES: Guide VI, section III.L.2 , Figure 6
