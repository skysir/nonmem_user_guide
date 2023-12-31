


 +--------------------------------------------------------------------+
 |                                                                    |
 |                   MTIME MNEXT MPAST MNOW MTDIFF                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Variables related to the model event time feature of PREDPP.
 CONTEXT:  Abbreviated code, verbatim code, user-supplied routines, NM-
 TRAN

 USAGE:
 MTIME(1)=THETA(1)

 Model event times are additional PK parameters defined in the PK  rou-
 tine or $PK block.  A model event time is not associated with any com-
 partment, but, like an absorption lag time, defines a  time  to  which
 the system is advanced.  When the time is reached, indicator variables
 are set and a call to PK is made.  At this call (and/or subsequent  to
 this  call)  PK or DES or AES or ERROR can use the indicator variables
 to change some aspect of the system, e.g., a term  in  a  differential
 equation,  or  the rate of an infusion.  This feature may be used with
 any ADVAN routine.  Model times  are  independent  of  non-event  dose
 times.   The following are reserved variables when used in abbreviated
 code.

 MTIME(i)
      The i-th model event time.  The maximum  number  of  model  event
      times is given by constant PCT in file SIZES (See sizes).

      MTIME(i)  may be less than, equal to, or greater than MTIME(i+1).
      Any MTIME(i) may be negative or have the value 0 (in  which  case
      MPAST(i)=1  always  and  MNEXT(i)=0  always;  see  below).  If PK
      defines MTIME(i) and MTIME(i+2) but not MTIME(i+1), then this has
      the  same  effect  as  defining  MTIME(i+1)=0.   PK  may redefine
      MTIME(i).  An ETA may be used in the definition of MTIME(i).

 MTDIFF
      The value of MTDIFF is 0 when PK is called.  If PK sets MTDIFF to
      a value other than 0, e.g., MTDIFF=1, then PREDPP will understand
      that with that call to PK, the values  of  one  or  more  of  the
      MTIME(i)  have  possibly  been reset.  It is not necessary to set
      MTDIFF at a call to PK with the first record of an individual  or
      with a reset record.

 The  following are the read-only indicator variables.  They are not to
 be set by PK.

 MNOW MNOW=i if MNEXT(i)=1 for some i.  MNOW=0 otherwise.

 MNEXT(i)
      MNEXT(i)=1 during the advance from the previous time to MTIME(i).
      Otherwise, MNEXT(i)=0.  The previous time may be an event time, a
      non-event time, or a model event time.

 MPAST(i)
      MPAST(i)=0 until the call to PK subsequent to the one  for  which
      MNEXT(i)=1.  At that call MPAST(i)=1. It then retains this value,
      unless MTIME(i) is redefined, in which case MPAST will be  appro-
      priately redefined as another step function.

 (See model time examples)
 (See Model Event Time: MNOW,MPAST,MNEXT)
 (See Model Event Time: MTIME)
 (See Model Event Time: MTDIFF)
 (See Circadian Example: Examples Using MTIME to Model Periodic Discon- |
 tinuities in $DES)

 REFERENCES: Guide VI, section III.F.9 
