


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               MODEL                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: MODEL subroutine
 CONTEXT: User-supplied subroutine; for use with PREDPP

 USAGE:
 SUBROUTINE MODEL(IDNO,NCM,NPAR,IR,IATT,LINK)
 USE PRMOD_CHAR, ONLY: NAME
 USE SIZES,     ONLY: DPSIZE,ISIZE
 INTEGER(KIND=ISIZE) :: IDNO,NCM,NPAR,IR,IATT,LINK
 DIMENSION :: IATT(IR,*),LINK(IR,*)

 GLOBAL DECLARATION:

 Versions before NONMEM 7.4:

 CHARACTER(LEN=8) :: NAME(PC)

 Versions after NONMEM 7.4:

 USE SIZES,     ONLY: SD
 CHARACTER(LEN=SD) :: NAME(PC)

 DISCUSSION:
 The  MODEL  subroutine is called by PREDPP only once at the start of a
 run when a general ADVAN (ADVAN 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18)
 is  used.   It  allows  the  user to specify aspects of the particular
 model he wishes to use.  When NM-TRAN is used, the $MODEL record  sup-
 plies this information.

 Input argument: None.

 Output argument:

  IDNO
      An  identification  number  for  the  MODEL  routine.   The value
      assigned by MODEL to IDNO is printed on the first PREDPP  problem
      summary page.

  NCM
      The total number of compartments in the model, excluding the out-
      put compartment.  Contains 0 when MODEL is called.  NCM  must  be
      no greater than PC-1.

  NPAR
      The  number  of basic PK parameters used in the PK routine.  Con-
      tains 0 when MODEL is called.  NPAR must be no greater than  con-
      stant  PG  in  file  SIZES.   With  the general non-linear models
      (ADVAN6, ADVAN8,  ADVAN9,  ADVAN13,  ADVAN14,  ADVAN15,  ADVAN16,
      ADVAN17, ADVAN18), NPAR may remain 0.

  IATT
      Values  of compartment attributes. The values of IATT(I,*) refers
      to the ith compartment.
      IATT(I,1)=0  initially off
      IATT(I,1)=1  initially on
      IATT(I,2)=0  may not be turned on and off
      IATT(I,2)=1  may be turned on and off
      IATT(I,3)=0  may not receive doses
      IATT(I,3)=1  may receive doses
      IATT(I,4)=0  not the default observation compartment
      IATT(I,4)=1  the default observation compartment
      IATT(I,5)=0  not the default dose compartment
      IATT(I,5)=1  the default dose compartment

      The remainder is used only with ADVAN9, ADVAN15, and ADVAN17:

      IATT(I,8)=0  not an equilibrium compartment
      IATT(I,8)=1  an equilibrium compartment
      IATT(I,9)=0  should not be included in the total drug amount
      in the system interior
      IATT(I,9)=1  should be included in the total drug amount
      in the system interior

  NAME
      Labels for the compartments (to be printed on the PREDPP  summary
      page under "FUNCTION").
      NAME(I)  is  the  label  for compartment i.  With versions before
      NONMEM 7.4, NAME was limited to 8  characters.   With  NONMEM  74
      NAME may be up to SD characters long.  SD is defined in SIZES and
      is set to 30 with NONMEM 7.4.

  LINK
      Only used with a general linear model (ADVAN5 and ADVAN7).
      LINK(I, J)=0: no drug may distribute from compartment i  to  com-
      partment j.
      LINK(I,  J)=K: drug distributes from compartment i to compartment
      j.  The rate constant quantifying this first  order  distribution
      is computed by PK and stored in the kth row of GG.
      When MODEL is called, all elements of LINK are 0.

 (See sizes).

 An initial steady state may be requested by routine MODEL;
 (See Initial Steady State: I_SS,ISSMOD)
 (See i_ss, initial_condition).

 REFERENCES: Guide IV, section V.C.4 
 REFERENCES: Guide VI, section VI.B , Figures 28, 29, 43
 REFERENCES: Guide VI, Appendix II 
