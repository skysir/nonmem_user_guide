


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          PASTZERO (NM75)                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE PRDDESLVU,    ONLY: PASTZERO

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL (KIND=DPSIZE) PASTZERO

 DISCUSSION:

 This  variable  is reserved with ADVAN18, for use with Delay Differen-
 tial Equations.  Sometimes you may wish to have the equations  transi-
 tion  from  the  past  to  the present other than at time T=0. In this
 case, set PASTZERO to a  non-zero  (including  negative)  value.   For
 example,

 $DES PASTZERO=-10.0

 REFERENCES: Guide Introduction_7
