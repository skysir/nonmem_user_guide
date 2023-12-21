


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           COMACT,COMSAV                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY: COMACT,COMSAV

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: COMACT,COMSAV

 DISCUSSION:

 COMACT is set by NONMEM; COMSAV is set by PRED.

 COMACT
      COMACT=0:  NONMEM is not making a copying pass, i.e., values from
      NMPRD4 will not be copied for tables and  scatterplots.   (NONMEM
      only  makes  a copying pass when PRED-defined items are listed in
      $TABLE or $SCATTER records.)
      COMACT=1: This is a copying pass  with  final  thetas;  etas  are
      zero.
      COMACT=2:  This  is  a  copying pass with final thetas and condi-
      tional eta estimates.
      COMACT=3: This is a copying pass with conditional (nonparametric)
      estimates  of  etas.   Such  a  pass takes place when the control
      stream includes this record:
       $NONPARAMETRIC ETAS

 COMSAV
      Set by PRED at ICALL<=1 or at COMACT=1 with the very  first  data
      record.  The Save Region is an initial part of NMPRD4.  If with a
      given data record, the value of a PRED-defined variable is stored
      in  this  region,  then  during  a copy pass with COMACT = 2, the
      value of the variable is initialized to the value obtained during
      the  copy  pass with the same data record and COMACT=1. COMSAV is
      the size of this region, i.e., the number of variables whose val-
      ues are stored in this region.

 When NM-TRAN is used, COMACT may be tested in abbreviated code. COMSAV
 may not be referenced in abbreviated code.  Instead, the COMSAV option
 of  the $ABBREVIATED record is coded, and NM-TRAN causes the specified
 value to be stored in COMSAV.

 Location prior to NONMEM 7: nmprd3

 REFERENCES: Guide IV, section III.B.7 , IV.E.2 
 REFERENCES: Guide VI, section III.J , IV.E 
