


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                DES                                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: DES subroutine
 CONTEXT:   User-supplied  subroutine;  for  use  with  PREDPP's  ADVAN
 6,8,9,13,14,15,16,17,18

 USAGE:

 SUBROUTINE DES (A,P,T,DADT,IR,DA,DP,DT)
 USE SIZES, ONLY: DPSIZE,ISIZE
 INTEGER(KIND=ISIZE) :: IR
 REAL(KIND=DPSIZE) :: A,P,DADT,DA,DP,DT
 DIMENSION :: A(*),P(*),DADT(*),DA(IR,*),DP(IR,*),DT(*)

 DISCUSSION:
 The DES subroutine is called by PREDPP to evaluate right-sides of dif-
 ferential equations.

 Input argument:

 A(n) The amount in the nth compartment at time T.

 P(n) The value of the nth PK parameter.

 T    Time. T takes values continuously over an integration interval.

 Output argument:

 DADT(n)
      The derivative with respect to T of the nth compartment's amount.
      It is important to note that PREDPP itself adds in the rates  for
      any infusions that may be active.

      It  is  possible to introduce drug into a compartment by explicit
      terms in a differential equation,  rather  than  by  PREDPP  dose
      event records.  Drug introduced in this manner is not included by
      PREDPP in the computation of the  output  compartment.   Specifi-
      cally,  the amount in the output compartment may be thought of as
      being calculated by summing all relevant  doses  from  the  INPUT
      file  (i.e.,  those  that precede the time of the present record,
      accounting for bioavailability), subtracting all amounts  present
      in  compartments other than the output compartment, and then mul-
      tiplying the result by the output fraction parameter.

      For example, suppose differential equations were used for ADVAN2,
      rather than the analytic solution.  They would be:
      DADT(1)=-P(3)*A(1)
      DADT(2)= P(3)*A(1)-P(1)*A(2)

 DA(n,j)
      The  derivative  of  DADT(n) with respect to A(j). Continuing the
      above example,
      DA(1,1)=-P(3)
      DA(2,1)= P(3)
      DA(2,2)=-P(1)

 DP(n,j)
      The derivative of DADT(n) with respect to  P(j).  Continuing  the
      above example,
      DP(1,3)=-A(1)
      DP(2,3)= A(1)
      DP(2,1)=-A(2)

 DT(n)
      The derivative of DADT(n) with respect to T.

 The  format  of arrays DA, DP, DT described above is called "full for-
 mat".  Alternately, compact format may be used  (See prdde1).   It  is
 the  default format for these arrays.  A description of compact format
 is beyond the scope of this document, but is described in  the  NONMEM
 7.4 version of Guide VI, Appendix IV. Compact Arrays in DES.

 Also  see variables in NONMEM modules, NONMEM-PRED modules, and PREDPP
 modules.
 (See variables in modules)
 In particular,
 (See DES AES: ICALL,IDEFD,IDEFA)

 REFERENCES: Guide IV, section V.C.7 
 REFERENCES: Guide VI, section VI.C 
