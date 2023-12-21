


 +--------------------------------------------------------------------+
 |                                                                    |
 |                       PRED-DEFINED VARIABLES                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE NMPRD4, ONLY: COM=>VRBL

 GLOBAL DECLARATION:
      MODULE NMPRD4
      USE SIZES, ONLY: DPSIZE
      IMPLICIT NONE
      SAVE
       REAL(KIND=DPSIZE), ALLOCATABLE, TARGET :: VRBL(:)
      END MODULE

 DISCUSSION:

 NMPRD4  is  a  module for PRED-defined variables (including PK-defined
 and ERROR-defined variables).  It contains an  array  VRBL.   (COM  is
 suggested above as an alternate to the name VRBL to be consistent with
 generated code, but this is optional.)  The VRBL  array  is  allocated
 dynamically  by NONMEM.  The default size is given by constant LNP4 in
 resource/sizes.f90, but this can be changed using the $SIZES record.

 Variables are stored in NMPRD4 for communication with other blocks  of
 abbreviated  code  or  with  user-written codes, or so that NONMEM has
 access to these variables  for  display  in  tables  or  scatterplots.
 (See $table, $scatter).   Values are stored by PRED and by subroutines
 of PREDPP. A SAVE region of NMPRD4 may be designated; (See COMACT,COM-
 SAV).

 When  abbreviated  code  is  used,  by default NM-TRAN lists all PRED-
 defined variables in the module, as well  as  certain  NM-TRAN-defined
 variables which give the values of partial derivatives.  Comment lines
 in the generated code (file FSUBS)  identify  these  variables.   This
 makes the listed variables globally accessible to all blocks of abbre-
 viated code, which may not be desirable.  (An exception  is  made  for
 initialization-finalization  variables;  these  are PRED-defined vari-
 ables that are first defined within a block that tests for ICALL 0, 1,
 3.   Such variables are stored locally and will not be changed by NON-
 MEM at ICALL 2 or 4.)  (See INFN-defined_variables.)

 Variables from a given block of abbreviated code can be excluded  from
 the module by including the following pseudo-statement in that block:
 COMRES=-1
 Variables from all blocks of abbreviated code can be excluded from the
 module by including the identical option (COMRES=-1) on the  $ABBREVI-
 ATED record.

 Variables  from  NMPRD4  can be used in a user-written subroutine.  If
 the user-written subroutine contains the USE NMPRD4 statement as shown
 above,  the  variables  have  the  names  COM(1),  COM(2), etc.  E.g.,
 X=COM(1)/COM(1).  It is possible to give the variables the same  (non-
 subscripted)  names  that  they  had  in $PRED, $PK or $ERROR. Suppose
 variables CL and V are used in $PK, and the following  is  present  in
 the generated PK routine:
      CL=>COM(00001); V=>COM(00002)

 In the user-written subroutine, after the last declaration, include:
      POINTER :: CL,V

 Prior to the first executable statement, include:
      CL=>COM(1); V=>COM(2)

 Now the code is X=CL/V.

 Code  generated  by NM-TRAN is somewhat more complex, but achieves the
 same result.

 REFERENCES: Guide III, section V.3.3 
 REFERENCES: Guide IV, section  III.B.7  ,  III.B.16  , III.B.17 , IV.H 
, V.C.5-9 
 REFERENCES: Guide VI, section III.J , III.L.2 , IV.E , IV.G.4  Figures 
5, 14
