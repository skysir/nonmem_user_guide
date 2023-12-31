


 +--------------------------------------------------------------------+
 |                                                                    |
 |                           CCONTR: F,G,H                            |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 NONMEM  code and previous examples that may be available from advanced
 users.

 USAGE:
      USE ROCM_REAL, ONLY: F=>FL2,G=>GL2,G2=>GGL2,H=>HL2

 GLOBAL DECLARATION:
      USE SIZES, ONLY: NO,LVR,LVR2,DPSIZE
      REAL(KIND=DPSIZE) :: FL2(NO),GL2(NO,LVR),GGL2(NO,LVR2*(LVR2+1)/2), &
                           HL2(NO,LVR*LVR/4+LVR/2)

 DISCUSSION:

 These variables change values with each L2 record.  They may  be  used
 by CCONTR.

  F   F(n)  =  value returned in F from PRED for the nth observation of
      the L2 record.

  G   G(n,i) = Partial derivative of F(n) with respect to eta(i)

  G2  G2(n,i*(i-1)/2+j)  =  second  partial  derivative  of  F(n)  with
      respect to eta(i) and eta(j) (j<=i)
      G2 is arranged in symmetric storage, e.g.

      1
      2   3
      4   5   6

      G2(n,1) = 2nd. partial of F(n) wrt. eta(1) eta(1)
      G2(n,2) = 2nd. partial of F(n) wrt. eta(2) eta(1)
      G2(n,3) = 2nd. partial of F(n) wrt. eta(2) eta(2)
      G2(n,4) = 2nd. partial of F(n) wrt. eta(3) eta(1)
      G2(n,5) = 2nd. partial of F(n) wrt. eta(3) eta(2)
      G2(n,6) = 2nd. partial of F(n) wrt. eta(3) eta(3)

  H   H(n,i)  =  Partial  derivative of F(n) with respect to eps(i) for
      i=1,neps, where neps is the number of epsilons in the problem.
      H(n,j*neps+i) = second partial derivative of F(n) with respect to
      eps(i) and eta(j).

      E.g. suppose that there are two epsilons in the problem:

      1   3   5
      2   4   6

      H(n,1) = partial derivative of F(n) wrt. eps(1)
      H(n,2) = partial derivative of F(n) wrt. eps(2)
      H(n,3) = 2nd. partial of F(n) wrt. eps(1) eta(1)
      H(n,4) = 2nd. partial of F(n) wrt. eps(2) eta(1)
      H(n,5) = 2nd. partial of F(n) wrt. eps(1) eta(2)
      H(n,6) = 2nd. partial of F(n) wrt. eps(2) eta(2)

 Location prior to NONMEM  7: rocm5

 REFERENCES: None.
