


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      PRIOR SIMULATION: ICMAX                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRIOR routine

 USAGE:
      USE NMPR_INT, ONLY: ICMAX=>IMAXSIM

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IMAXSIM

 DISCUSSION:

 This variable allows the user to control the behavior of the NWPRI and
 TNPRI NONMEM utility routines.  It is relevant during simulation  when
 PRIOR  sets  ISS  to  a  non-zero  value.  By default, NWPRI and TNPRI
 attempt at most 100 times to obtain a sample.  After  that  number  of
 attempts, they terminate with the message

 MAXIMUM ATTEMPTS TO OBTAIN SAMPLE EXCEEDED

  ICMAX
      If  a  value  is set in ICMAX, this overrides the default of 100.
      E.g., in a PRIOR subroutine:

       USE NMPR_INT, ONLY: ICMAX=>IMAXSIM
   ... other code ...
       ICMAX=500

 Location prior to NONMEM 7: nmpr50

 REFERENCES: None.
