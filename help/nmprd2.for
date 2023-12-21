


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         PRED ERROR MESSAGE                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_CHAR, ONLY :: ETEXT

 GLOBAL DECLARATION:
      CHARACTER(LEN=132) :: ETEXT(3)

 DISCUSSION:

 ETEXT
      The user message which appears as part of the PRED error message.
      When PRED returns to NONMEM with IERPRD>0, PRED  may  store  here
      one  to  three  lines of text explaining the error, for printing.
      The number of lines of text to  be  printed  must  be  stored  in
      NETEXT. (See PRED Exit Code)

 When  abbreviated  code is used, the EXIT statement causes appropriate
 text to be stored in ETEXT.

 Location prior to NONMEM 7: nmprd2

 REFERENCES: Guide VI, section III.K.2 , IV.F 
 REFERENCES: Guide IV, section IV.G 
