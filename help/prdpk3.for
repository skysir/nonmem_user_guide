


 +--------------------------------------------------------------------+
 |                                                                    |
 |                  COMPARTMENT INITIALIZATION: A_0                   |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: PK routine

 USAGE:
      USE PRMOD_REAL, ONLY: A_0,DA_0,D2A_0

 GLOBAL DECLARATION:
      USE SIZES, ONLY: PC,LVR,DPSIZE
      REAL(KIND=DPSIZE) :: A_0(PC),DA_0(PC,LVR),D2A_0(PC,LVR,LVR)

 DISCUSSION:

 PREDPP  sets A_0FLG to 1 (See compartment initialization: a_0flg) at a
 call to PK with the first event record of an individual record (if the
 data are population data), with the first event record of the data set
 (if the data are single-subject data), and with a  reset  record.   At
 such  times, the amounts in the various compartments can be set by the
 PK routine.  It can do this by storing  the  initial  values  for  the
 state  vector  and its partials in A_0,DA_0,D2A_0.   The amount in the
 output compartment can not be set.

 A_0
      A_0(n) = the amount for compartment n

 DA_0
      DA_0(n,i) = the derivative of A_0(n) wrt eta(i)

 D2A_0
      D2A_0(n,i,j) = the second derivative of A_0(n) wrt eta(i), eta(j)
      (lower-triangular; j=1, ..., i)

 There  is  a  one-to-one  correspondence  between  A_0,DA_0,D2A_0  and
 A,DAETA,D2AETA (See State Vector: A).

 NM-TRAN includes A_0,DA_0,D2A_0 in the PK routine when the  $PK  block
 includes  references  to  variables A_0FLG, A_0, or A_INITIAL, or when
 verbatim code is present.

 (See Compartment Initialization: A_0FLG, State Vector: A).

 Location prior to NONMEM 7: prdpk3

 REFERENCES: None.
