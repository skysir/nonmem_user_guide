


 +--------------------------------------------------------------------+
 |                                                                    |
 |           NPD,NPDE,NPDE_MODE,DV_LOQ,CDF_L,DV_LAQ_CDF_LA            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED routine

 USAGE:
      USE NMPRD_INT, ONLY: NPDE_MODE
      USE NMPRD_REAL, ONLY: DV_LOQ

 GLOBAL DECLARATION:
             INTEGER(KIND=ISIZE) ::  NPDE_MODE
             REAL(KIND=DPSIZE) :: DV_LOQ

 DISCUSSION:

 NONMEM 7 can compute the following:

 NPDE
      Monte-Carlo generated normalized probability distribution error.

 NPD  The correlated (or non-decorrelated) NPDE value.

 These are calculated during the Table Step, when NPD or NPDE is listed
 in $TABLE or $SCATTER.

 NONMEM  sets  NPDE_MODE=1  when  it  is  computing  them.   Otherwise,
 NPDE_MODE=0.

 EXAMPLE:

 To  incorporate  LOQ  (Limit of Quantification) data into NPDE assess-
 ments [4], use the following method (as an example).
 See INTRODUCTION TO NONMEM 7, Reference [4].
 Here, TYPE and LOQ are user-defined in previous code, or data item.

 $ERROR
 SD = THETA(5)
 IPRED = LOG(F)
 DUM = (LOQ - IPRED) / SD
 CUMD = PHI(DUM)
 IF (TYPE .EQ. 1.OR.NPDE_MODE.EQ.1) THEN
       F_FLAG = 0
       Y = IPRED + SD * ERR(1)
 ENDIF
 IF (TYPE .EQ. 2.AND.NPDE_MODE.EQ.0) THEN
       F_FLAG = 1
       Y = CUMD
       MDVRES=1
 ENDIF
 IF(TYPE.EQ.2) DV_LOQ=LOQ
   ....
 $TABLE NPD NPDE

 By default, DV_LOQ is set to -1.0d-300  by  the  NONMEM  routine  that
 calls PREDPP/PRED.  If the user's ERROR/PRED sets DV_LOQ to some other
 value and NPDE_MODE=1, then the NPDE is being  evaluated  during  that
 time,  and this censored value is to be treated as if it is a non-cen-
 sored datum with value of LOQ (DV_LOQ=LOQ), in  accordance  with  [4],
 utilizing  a  standard  F_FLAG=0  definition  for Y.  Note that during
 estimation of the objective function (when NPDE_MODE=0), NPDE  is  not
 being evaluated, and censored values should be treated using F_FLAG=1,
 and Y must be defined as the integral of the normal density from  -inf
 to LOQ.

 New in nm74, for use with NPD, the user may supply the cumulative dis-
 tribution function using the reserved variable CDF_L.

 New in nm743, you can specify an above  quantifiable  limit  with  the
 reserved  parameter DV_LAQ as well.  The CDF reserved variable associ-
 ated with above quantitation level DV_LAQ is CDF_LA.

 The DV_LOQ/CDF_L and DV_LAQ/CDF_LA reserved variables may be also used
 for evaluating NPD diagnostics for categorical data.

 See INTRODUCTION TO NONMEM 7, MDVRES=0

 REFERENCES: Guide Introduction_7
