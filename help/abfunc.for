


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        ABBREVIATED FUNCTION                        |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: User-supplied function that may be used in abbreviated code.
 CONTEXT: Fortran coded function

 The following applies to all versions of NONMEM.

 USAGE:
 FUNCTION FUNCA(X,X1,X2)
 DOUBLE PRECISION X,X1,X2,FUNCA
 DIMENSION X(9),X1(9),X2(9,9)

 The constant "9" must be used exactly as shown.

 DISCUSSION:
 The  FORTRAN  function  FUNCA  may be used in abbreviated code.  Simi-
 larly,
 the functions FUNCB, FUNCC, etc.  may  be  used;  their  constructions
 would  be similar to that of FUNCA.  In abbreviated code, the function
 is referenced with a single argument, either a (scalar) expression  or
 a vector.  Certain reserved names such as VECTRA, VECTRB, etc.  may be
 used in the abbreviated code, but the name of the argument is  not  in
 the  code defining the function.  (See abbreviated code).  In the code
 defining the function, the function has three arguments.

 Input Argument:

  X   The value of the argument (which may be a vector).

 Output Argument:

  FUNCA
      The value of the function is to be stored in FUNCA.

 There are two other outputs.

  X1  X1(n) is  the  first-partial  derivative  of  the  function  with
      respect to the nth element of the argument.
      If  the  value  of X(n) will not be a value of a random variable,
      X1(n) need not be set.
      If the argument X is a (scalar) expression, only X1(1) is needed.

  X2  X2(n,m) is the second-partial derivative  of  the  function  with
      respect to the nth and mth elements of the argument.
      If the value of either X(n) or X(m) will not be a value of a ran-
      dom variable, X2(n,m) need not be set.  Nor need any value of  X2
      be set if a test of MSEC=1 is false (See Partial Derivative Indi-
      cators).
      If the argument is a (scalar) expression,  only  X2(1,1)  may  be
      needed.   If the argument is a vector and if the value of X2(n,m)
      is needed, then the value of X2(m,n)  is  needed  as  well,  even
      though these two values will be identical.

 The  Fortran  code  for  the functions should be placed in one or more
 files.  Suppose there is one such file and its name  is funcfile.   It
 should be listed on the $SUBROUTINES record.  It may contain more than
 one FUNCTION.  E.g.,
 $SUBROUTINES ... OTHER=funcfile

 With NONMEM VI 2.0, reserved function names are FUNCA  through  FUNCC.
 reserved names for vectors are VECTRA through VECTRC.

 With NONMEM 7.3, reserved function names are FUNCA through FUNCI.

 With  NONMEM  7.4,  reserved  function  names are FUNCA through FUNCZ.
 Reserved names for vectors are VECTRA through VECTRZ.

 If the abbreviated code uses one of these reserved function names, but
 the  user  does  not  supply the code via $SUBR OTHER, the NONMEM exe-
 cutable will be created, but there will be an  error  message  in  the
 NONMEM report file such as:
 FUNCA WAS CALLED, BUT NO CODE WAS SUPPLIED.

 With  NONMEM  7.4,  extended  reserved names are recognized. These are
 FUNCxy and FUNCxyz, where each of x, y, z  stands  for  an  alphabetic
 character  A-Z,  e.g.,  FUNCAB  or FUNCABC.  Similar extended reserved
 names for vectors are also recognized: e.g, VECTRAB or VECTRABC.
 If the abbreviated code uses one of these extended  reserved  function
 names, but the user does not supply the code via $SUBR OTHER, the NON-
 MEM executable cannot be created, and there will be an  error  message
 from the Fortran compiler such as:
 Undefined Symbols
 _funcaaa_

 (See Abbreviated function example).

 The following applies only to NONMEM 7.4 and higher.

 The  $ABBR  FUNCTION  allows the user to declare the name of the func-
 tion, the name of the vector of input arguments (which is of no use to
 the  function), and the length of the vector of input arguments, which
 is passed to the function as argument NDIM.  Constant NVECX  in  SIZES
 gives the maximum value size of the vector of input arguments.

 In the NONMEM control stream:

 $ABBR FUNCTION function_name(input_vector_name,dimension,usage)

 (See $ABBREVIATED).
 (See VECTORS and ABBREVIATED FUNCTIONS)

 The function should be coded as follows:

 FUNCTION function_name(X,X1,X2,NDIM)
 USE SIZES, ONLY: DPSIZE,ISIZE
 INTEGER(KIND=ISIZE) :: NDIM
 REAL(KIND=DPSIZE), INTENT(IN)     :: X
 REAL(KIND=DPSIZE), INTENT(IN OUT) :: X1,X2
 REAL(KIND=DPSIZE):: function_name
 DIMENSION X(NDIM),X1(NDIM),X2(NDIM,NDIM)

 Function_name  must  be  exactly as specified via $ABBR FUNCTION.  The
 value coded as "dimension" is passed to the function  as  NDIM.   (The
 value coded as "usage" is optional and does not appear in code for the
 function.)

 The value of the function is to be stored  in  function_name.   Other-
 wise, the function is coded exactly as described above, with X1 and X2
 as output arguments. The input_vector_name is never used in  the  For-
 tran code for the function.

 REFERENCES: none
