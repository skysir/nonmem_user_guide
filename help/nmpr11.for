


 +--------------------------------------------------------------------+
 |                                                                    |
 |                     SIMULATION: IETAOL IEPSOL                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine
 USAGE:
      USE NMPR_INT, ONLY: IETAOL,IEPSOL

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: IETAOL,IEPSOL

 DISCUSSION:

 Values  are stored by PRED for use with SIMETA and SIMEPS.  They allow
 the effect of the NEW option on the $SIMULATION record to be  overrid-
 den,  which  may  be  useful when the PRED repetition feature is used.
 They also allow calls to SIMETA and SIMEPS to be ignored.  This may be
 useful  when  e.g.  a  simulation  was undertaken with user code where
 SIMEPS is called with every data record (as happens automatically with
 NM-TRAN  generated  codes), and the exact same simulation is now to be
 repeated, but with a data set obtained from the  earlier  one  by  the
 addition  of  new  nonobservation records (with which SIMEPS output is
 not needed).  If calls to SIMEPS with the new records are not ignored,
 SIMEPS  output will be generated with all the records, and in particu-
 lar, SIMEPS output will be generated with the original records,  which
 will differ from what it was earlier, thus resulting in the simulation
 of a different set of observations.

 IETAOL

      IETAOL=-1: Next call to SIMETA will be ignored

      IETAOL=0: Next call to SIMETA will behave as usual.

      IETAOL=1: Next call to SIMETA will behave as though NEW  had  not
      been  specified.   If  SIMETA has been called previously with the
      individual record, SIMETA will produce the previous eta values.

 IEPSOL

      IEPSOL=-1: Next call to SIMEPS will be ignored.

      IEPSOL=0: Next call to SIMEPS will behave as usual.

      IEPSOL=1: Next call to SIMEPS will behave as though NEW  had  not
      been  specified.   If  SIMEPS has been called previously with the
      level-two record, SIMEPS will produce the previous  epsilon  val-
      ues.

 Values must be stored before PRED calls SIMETA (SIMEPS).  NM-TRAN gen-
 erated or Library code has a call to SIMETA  (SIMEPS)  in  its  second
 section  (see  Guide IV).  If this call is to be affected, values must
 be stored using verbatim code in the FIRST block.

 Location prior to NONMEM 7: nmpr11.for

 REFERENCES: None.
