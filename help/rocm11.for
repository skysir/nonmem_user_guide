


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    MIXTURE MODEL: MIXNUM,MIXEST                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_INT,  ONLY: MIXNUM=>MIXCALL,MIXEST=>IMIXEST

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: MIXCALL,IMIXEST

 DISCUSSION:

  MIXNUM
      The number of the subpopulation for which PRED is to compute out-
      puts.

  MIXEST
      The number of the subpopulation estimated to be that  from  which
      the  individual's  data  most probably arises.  Should be used at
      ICALL = 3 and at ICALL = 2 when COMACT > 0.

 Used with mixture models.  MIXNUM and MIXEST are reserved variables in
 all blocks of abbreviated code.

 This is a general description. For details, see mixnum.                |

 (See MIXNUM_MIXEST_MIXP)
 (See Mixture_model:_MIXP, mix)
 (See mixture_model_example)

 Location prior to NONMEM 7: rocm11

 REFERENCES: Guide VI, section III.L.2 , Figure 6
