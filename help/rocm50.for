


 +--------------------------------------------------------------------+
 |                                                                    |
 |                OBJECTIVE FUNCTION VALUE INDIVIDUAL                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_INT, ONLY: IIDX=>IDVALX
      USE ROCM_REAL, ONLY: CNTID=>OFV_IND

 GLOBAL DECLARATION:
      USE SIZES, ONLY: MAXIDS,DPSIZE
      INTEGER(KIND=ISIZE) :: IDVALX(MAXIDS)
      REAL(KIND=DPSIZE) :: OFV_IND(MAXIDS)

 DISCUSSION:

 Note:  With NONMEM 7, the additional output file root.phi contains the
 same information.
 (See additional_output_file).

 These variables contain values of the ID data item and individual con-
 tributions  to  the  objective  function.   The values are in data-set
 order.

 IIDX Values of the ID data item.

 CNTID
      Values of the individual contribution to the  objective  function
      for the corresponding values of IIDX.

 E.g.,  IIDX(n) is the ID data item for the nth. individual record, and
 CNTID(n) is the contribution to the objective function  for  the  nth.
 individual record.

 These  values  should  only  be  displayed  at ICALL = 3 (finalization
 block).

 With NONMEM VI 1.0, they can only be displayed using verbatim code.
 With NONMEM VI 2.0 and later releases, they can be used on  the  right
 and  displayed  using abbreviated code in $PRED, $PK, $ERROR and $INFN
 blocks (See Individual objective function example).

 (See write).

 They may also be displayed in a table, using

 $ABBR COMRES=2

 and code such as the following in the $ERROR or $PK block:

 IF (COMACT.EQ.1) THEN
  COM(1)=IIDX(NIREC)                                                    |
  COM(2)=CNTID(NIREC)                                                   |
 ENDIF
 The following, for example,  will produced a separate  table  for  the
 values:

 $TABLE IID=COM(1) CNT=COM(2) FILE=comvals NOAPPEND NOPRINT FIRSTONLY

 Note: With earlier versions than NONMEM 7.3, verbatim code is needed:  |
 "  COM(1)=IIDX(NIREC)                                                  |
 "  COM(2)=CNTID(NIREC)                                                 |

 Location prior to NONMEM 7: rocm50

 REFERENCES: None.
