


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          STANDARD ERRORS                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: PRED routine

 USAGE:
      USE ROCM_REAL, ONLY: SETHET=>SETH,SEOMEG=>SEOM,SESIGM=>SESIG,
          SETHETR=>SETHR

 GLOBAL DECLARATION:
      USE SIZES, ONLY: LTH,LVR,DPSIZE
      REAL(KIND=DPSIZES) :: SETH(LTH),SEOM(LVR,LVR),SESIG(LVR,LVR)

 DISCUSSION:

 SETHET
      SETHET(i) = the standard error of the estimate of theta(i).

 SETHETR
      SETHETR(i)  =  the  standard  error  of  the estimate of reported |
      theta(i).  If record $THETAR is not used, SETHET and SETHETR  are |
      equal.   If  record  $THETAR is used, then SETHET is the standard |
      error of the internal theta as used in $PK/$PRED, and SETHETR  is |
      the standard error of the theta reported in the report file.

 SEOMEG
      SEOMEG(i,j) = the standard error of the estimate of omega(i,j).

 SESIGM
      SESIGM(i,j) = the standard error of the estimate of sigma(i,j).

 These values should only be used with ICALL = 3.

 Location prior to NONMEM 7: rocm7

 REFERENCES: Guide I, section C.3.5.2 
 REFERENCES: Guide V, section 5.4.2.1 
