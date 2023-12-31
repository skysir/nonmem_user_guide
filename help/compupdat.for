


 +--------------------------------------------------------------------+
 |                                                                    |
 |                  COMPARTMENT UPDATE BLOCK (NM75)                   |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for compartment update
 CONTEXT: $PK abbreviated code

 SAMPLE:
 $PK
 IF (A_UFLG.EQ.1) THEN
  ... compartment update block ...
 ENDIF

 DISCUSSION:
 With  NONMEM 7.5, a compartment update block is a block of abbreviated
 code that is very similar to a compartment initialization block.

 In a compartment initialization block, PREDPP sets A_0FLG to  1  at  a
 call  to  PK  with all the compartments at their initial state so that
 values may be assigned to reserved variables A_0(n).
 (See Compartment Initialization: A_0)

 In a compartment update block, the user sets A_UFLG  to  1  in  PK  to
 indicate  to  PREDPP that PK is going to update the compartments.  The
 desired compartment values may be set in the array A_U(n).   The  user
 should  use  MTIME  to  designate a variable time position at which an
 abrupt change in compartment amounts occurs.  One could input  a  dose
 as follows:

 MTIME(1)=wtime
 MTDIFF=1
 AZTEST=A_0FLG
 IF(TSTATE==MTIME(1).AND.AZTEST==0) A_UFLG=1
 IF(A_UFLG==1) THEN
 A_U(1)=A(1)+wdose
 A_U(2)=A(2)
 A_U(3)=A(3)
 ENDIF

 With  the Compartment Update Block, the user sets A_UFLG to 1 when the
 compartments are to be updated.  The A_UFLG event  must  be  triggered
 with  an IF(TSTATE==MTIME()) condition as indicated in the above exam-
 ple.  Values may be assigned to reserved variables A_U(n).  The  value
 of  the  amount  in the  nth compartment (the nth element of the state
 vector) is set to the  value  assigned  to  A_U(n).   Any  A_U(x)  not
 explicitly  defined are set to 0.  An un-assigned A_U(k) should retain
 its value, A_u(k)=A_u(k).

 The code "IF(A_UFLG==1)...THEN...ENDIF" is optional,  as  NMTRAN  will
 insert  it if not present.  A_0FLG must be 0 whenever A_UFLG is set to
 1, as shown in the example above (..\examples_uflg.ctl).

 The rules for compartment update blocks are similar to those for  com-
 partent initialization blocks.

 PREDPP  expects to find the A_U values in the A_0 arrays.  NMTRAN con-
 verts A_U() in abbreviated code to A_0() during FSUBS  code  construc-
 tion.

 (See  Guide  Introduction_7,  "Updating Amounts in Compartments at any
 Time: The A_UFLG Flag (NM75)"

 REFERENCES: Guide Introduction_7
