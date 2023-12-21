


 +--------------------------------------------------------------------+
 |                                                                    |
 |                          PASTZERO (NM75)                           |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP-PK global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE PRRADAR5U, ONLY: PASTZERO

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL (KIND=DPSIZE) PASTZERO

 DISCUSSION:

 This variable is reserved with ADVAN16 and ADVAN17, for use with Delay
 Differential Equations.  Sometimes you may wish to have the  equations
 transition  from  the  past  to the present other than at time T=0. In
 this case, set PASTZERO to a non-zero (including negative) value.  For
 example,

 $DES PASTZERO=-10.0

 REFERENCES: Guide Introduction_7
