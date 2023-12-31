


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                AES                                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: AES subroutine
 CONTEXT:  User-supplied  subroutine;  for use with PREDPP's ADVAN9 and
 ADVAN15 and ADVAN17

 USAGE:
 SUBROUTINE AES (INIT,A,P,T,E,IR,DA,DP,DT)
 USE SIZES, ONLY: DPSIZE,ISIZE
 INTEGER (KIND=ISIZE) :: INIT,IR
 REAL(KIND=DPSIZE) :: A,P,E,DA,DP,DT
 DIMENSION :: A(*),P(*),E(*),DA(IR,*),DP(IR,*),DT(*)

 DISCUSSION:
 The AES subroutine is called by PREDPP to evaluate  algebraic  expres-
 sions for ADVAN9, ADVAN15, and ADVAN17.

 Input argument:

 P(n) The value of the nth PK parameter.

 T    Time. T takes values continuously over an integration interval.

 Input/output arguments:

 INIT When AES is called with INIT=1, this is an initial condition call
      at the start of an integration interval.  Approximate amounts  in
      each  equilibrium compartment n at time T must be stored in A(n).
      (The amounts in  the  non-equilibrium  compartments  are  already
      stored in the lower-numbered elements of A.)  DA, DP, DT need not
      be computed.

      If AES stores approximate values in A, it must set INIT=0.
      If AES stores exact values in A, it must leave INIT unchanged.

      When AES is called with  INIT=2,  the  values  of  the  algebraic
      expressions  must  be  stored in E and the derivatives in DA, DP,
      and DT.

 A(n) The amount in the nth compartment at time T.

 Output argument:

 E(k) The value of the algebraic expression g(k).

 DA(k,j)
      The derivative of g(k) with respect to A(j).

 DP(k,j)
      The derivative of g(k) with respect to P(j).

 DT(k)
      The derivative of g(k) with respect to T.

 Also see variables in NONMEM modules, NONMEM-PRED modules, and  PREDPP
 modules.
 (See variables in modules)
 In particular,
 (See DES AES: ICALL,IDEFD,IDEFA)

 REFERENCES: Guide IV, section V.C.8 , V.C.9 
 REFERENCES: Guide VI, section VI.E 
