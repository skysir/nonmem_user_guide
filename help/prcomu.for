


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      BIOAVAILABILITY BEHAVIOR                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP global variables
 CONTEXT: For use with PREDPP

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 PREDPP code.

 USAGE:
      $PK
      "FIRST
      " USE PRCOM_LOG, ONLY: NEWBIO=>NEWWAY
      "MAIN
      " NEWBIO=.FALSE.

 GLOBAL DECLARATION:
      LOGICAL NEWWAY

 DISCUSSION:
 This  variable allows the user to change the way bioavailability frac-
 tion parameters are used by PREDPP.  Prior to 1990 (i.e., with  NONMEM
 III  /  PREDPP  I) bio-availability fractions applied to all transient
 doses except infusions (doses with AMT > 0 and RATE > 0 ).  Similarly,
 it did not apply to steady-state doses with multiple infusion and with
 a positive rate data item.

 After the change, bio-availability fractions apply  to  all  transient
 (non  steady-state) doses, and to steady-state doses for which the AMT
 is specified.  The only doses to which bio-availability  fractions  do
 not apply are steady-state with constant infusion.

 By  default, NEWBIO is true and the new rules apply.  If NEWBIO is set
 to false, bioavailability fraction parameters revert to their pre-1990
 behavior.

 Location prior to NONMEM 7: prcomu

 REFERENCES: None.
