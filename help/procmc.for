


 +--------------------------------------------------------------------+
 |                                                                    |
 |                 COMPARTMENT INITIALIZATION: A_0FLG                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP read-only global variables
 CONTEXT: User-supplied PK routine

 USAGE:
      USE PROCM_INT, ONLY: A_0FLG

 GLOBAL DECLARATION:
      INTEGER(KIND=ISIZE) :: A_0FLG

 DISCUSSION:

 A_0FLG
      When  ICALL>=2,  at  a  call to PK with the first record of of an
      individual record or with a reset record, PREDPP sets  A_0FLAG=1.
      At  such  a call, the state vector A and its first and second eta
      partials (DAETA and D2AETA) have been set to zero.   A_0FLG=0  at
      all other calls to PK.

 When  A_0FLG=1,  PK  may  initialize  compartments.  It can do this by
 storing the initial values for the state vector and  its  partials  in
 A_0,DA_0,D2A_0.

 NM-TRAN  includes A_0FLG in the PK routine when the $PK block includes
 references to variables A_0FLG, A_0(n), or A_INITIAL(n), or when  ver-
 batim code is present.

 (See Compartment Initialization: A_0, State Vector: A).

 Location prior to NONMEM 7: procmc

 REFERENCES: None.
