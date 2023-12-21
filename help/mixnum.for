


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         MIXNUM_MIXEST_MIXP                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Variables used with mixture models
 CONTEXT: Abbreviated code and user-supplied subroutines

 USAGE:
      IF (MIXNUM.EQ.1) THEN
        ...
      ENDIF
      IF (MIXNUM.EQ.2) THEN
        ...
      ENDIF
      IF (MIXNUM.EQ.3) THEN
        ...
      ENDIF

 DISCUSSION:

 With mixture models, MIXNUM, MIXEST and MIXP are global variables that
 may be used as right-hand quantities (or  in  logical  conditions)  in
 abbreviated code or user-supplied routines.
 In general:
 MIXNUM is an input to PRED (or PREDPP) set by NONMEM.
 MIXEST is an output (result) or consequence from the estimation set by
 NONMEM.

 MIXNUM

 This is the index of the subpopulation for which variables are  to  be
 computed.

 During  ICALL=4,  MIXNUM  is  the index of the sub-population that was
 randomly selected to simulate the subject's data.   When  ICALL=3  and
 the PASS subroutine used, or when ICALL=5 or 6, MIXNUM is the index of
 the subpopulation estimated to be that  from  which  the  individual's
 data most probably arises, and is equal to MIXEST (see below).

 MIXEST

 When  ICALL=3,5,6,  and at ICALL=2 when COMACT is not 0, MIXEST is the |
 index of the subpopulation estimated to be that from which  the  indi- |
 vidual's data most probably arises.  Before nm7.3, at all other condi- |
 tions, the value of MIXEST is not meaningful, such as when ICALL=4  or |
 ICALL=2  and  COMACT=0.   As  of NONMEM 7.3, when ICALL=4, MIXEST will |
 equal MIXNUM, and when ICALL=2 and COMACT=0, MIXEST will be 0 to indi- |
 cate it has not been assigned a meaningful value.

 MIXP

 These are the mixture probabilities P(i) computed by subroutine MIX or
 by the $MIX block of abbreviated code.

 In abbreviated code, MIXP may be coded in either of three ways:

 MIXP(MIXNUM)
 MIXP
 MIXP(i)

 The first two ways are identical; i.e., when no  subscript  is  coded,
 NM-TRAN supplies the subscript MIXNUM.

 With  the third, the index i must not exceed the value of MMX in SIZES |
 (or as set by $SIZES record).

 (See Mixture_model:_MIXNUM,MIXEST)
 (See Mixture_model:_MIXP)
 (See mix, mixture model example).

 REFERENCES: Guide IV, section III.B.4 , III.B.6 
 REFERENCES: Guide IV, section IV.E.1 , 4.E.2 
 REFERENCES: Guide VI, section III.L.2 , Figure 6
