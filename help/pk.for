


 +--------------------------------------------------------------------+
 |                                                                    |
 |                                 PK                                 |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PK subroutine
 CONTEXT: User-supplied subroutine; for use with PREDPP

 USAGE:

 Versions before NONMEM 7.2:

 SUBROUTINE PK(ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,IRGG,GG,NETAS)
 USE SIZES,     ONLY: DPSIZE,ISIZE,LVR
 REAL(KIND=DPSIZE) :: EVTREC
 INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS,IRGG,NETAS
 DIMENSION :: IDEF(7,*),THETA(*),EVTREC(IREV,*),INDXS(*)
 DIMENSION :: GG(IRGG,LVR+1,*)

 With NONMEM 7.2 and higer:

 SUBROUTINE PK(ICALL,IDEF,THETA,IREV,EVTREC,NVNT,INDXS,IRGG,GG,NETAS)
 USE SIZES,     ONLY: DPSIZE,ISIZE
 USE PRDIMS,    ONLY: GPKD
 IMPLICIT REAL(KIND=DPSIZE) (A-Z)
 REAL(KIND=DPSIZE) :: EVTREC
 INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS,IRGG,NETAS
 DIMENSION :: IDEF(7,*),THETA(*),EVTREC(IREV,*),INDXS(*)
 DIMENSION :: GG(IRGG,GPKD+1,*)

 DISCUSSION:
 The  PK  subroutine is called by PREDPP to obtain values for basic and
 additional pharmacokinetic parameters.  Basic PK parameters are  typi-
 cally the rate constants ("micro-constants") for use in kinetic formu-
 las.  PK can compute instead parameters such as clearance and  volume,
 and  a translator ("TRANS") subroutine can be used to convert these to
 rate constants.  Additional PK parameters  include  compartment  scale
 parameters,  which  PREDPP uses to convert compartment amounts to con-
 centrations, and dose-related  parameters  such  as  modeled  infusion
 rates.

 Input argument:

  ICALL
      ICALL=1:  PK  has been called for initialization at the beginning
      of a NONMEM problem; one such call per problem.  EVTREC  contains
      the  first  event  record.  THETA contains the initial estimates.
      PK must return values in IDEF, which inform PREDPP what tasks  it
      will  perform  at  later calls.  It may also set GG (k,1,1) to 1,
      indicating that the kth PK parameter will be modeled as a log  (a
      feature which cannot be specified using abbreviated code).

      ICALL=2:  PK has been called to obtain parameter values; multiple
      calls occur.  PK should call GETETA to obtain ETA values for  the
      individual. PK should compute all relevant portions of GG.

      ICALL=4:  PK has been called during the Simulation Step; multiple
      calls occur.  PK should call SIMETA to obtain simulated ETA  val-
      ues  for  the  individual.  PK should compute individual-specific
      parameters for column 1 of GG.

      ICALL=5: PK has been called during the  computation  of  expecta-
      tions;  multiple calls occur.  Such calls occur when the marginal
      (MRG_) data item is defined in the data set and has non-zero val-
      ues  for  some records.  If the MRG_ data item has the value 1 or
      2, the values of PK-defined items contribute to the  expectations
      computed for these items.

      ICALL=6:  PK  has  been called during the computation of raw data
      averages; Such calls occur when the raw-data data item (RAW_)  is
      defined in the data set and has non-zero values for some records.

  THETA
      The NONMEM THETA vector.

  EVTREC
      The PREDPP event record.

  INDXS
      The  values specified in the $INDEX record of the NM-TRAN control
      stream.  (This is the NONMEM Index array starting at position 12,
      the first position beyond those positions used by PREDPP itself).

  NETAS
      The number of population etas in the problem.

 Output argument:

  IDEF
      PK should store values in IDEF only when ICALL=1.
      Values may be stored in the following elements:
      IDEF(1,1)=-9  (required)
      IDEF(1,2) is the call limiting element (compare $PK's CALLFL).
      Values are:
        -2:  call  with every event record and at additional and lagged
        dose times.
        -1: call with every event record (default)
         0: call with first event record and new TIME values
         1: call once per individual record

      Compartment initialization may be performed by routine PK;
      (See Compartment Initialization: A_0)
      (See sCompartment Initialization: A_0FLG).
      An initial steady state may be requested by routine PK;
      (See Initial Steady State: I_SS,ISSMOD)
      (See i_ss, initial_condition).
      The value in IDEF(1,3) describes whether PK performs  compartment
      initialization,  i.e.,  whether or not PK initializes elements of
      the initial state vector A_0(n).  Values are:
        -1: PK may initialize A_0.
         0: PK does not initialize A_0.
         1: PK does initialize A_0.
      The default used by PREDPP is IDEF(1,3)=-1.  However,  when  com-
      partment  initialization is not implemented, then if IDEF(1,3) is
      set to  0,  PREDPP  can  avoid  some  time-consuming  processing.
      Indeed,  when  $PK  abbreviated or verbatim code is supplied, and
      there is  no  reference  to  compartment  initialization  amounts
      A_0(n)  in  either the abbreviated or verbatim code, then NM-TRAN
      sets IDEF(1,3)=0.

      Compartment amounts may be used by routine PK;
      (See State Vector: A).
      The value in IDEF(1,4) describes whether PK uses  derivatives  of
      compartment amounts (e.g. compartment amounts themselves are used
      as random variables in arithmetic statements in PK).  Values are:
        -1: PK may use derivatives of compartment amounts.
         0: PK does not use derivatives of compartment amounts.
         1: PK uses derivatives of compartment amounts.
      The default used by PREDPP is IDEF(1,4)=-1.  However, when deriv-
      atives  of compartment amounts are not used, then if IDEF(1,4) is
      set to  0,  PREDPP  can  avoid  some  time-consuming  processing.
      Indeed,  when  $PK  abbreviated or verbatim code is supplied, and
      there is no reference to A(n) (as a random variable in an  arith-
      metic  statement)  in  the abbreviated code (or to derivatives of
      A(n) in the verbatim code), then NM-TRAN sets IDEF(1,4)=0.

      Remaining elements contain row numbers in GG:
      IDEF(2,1): row number of F0 (output fraction)
      IDEF(2,2): row number of XSCALE (Time Scale)
      IDEF(2,3): row number of lowest-numbered MTIME
      IDEF(2,4): row number of highest-numbered MTIME
      IDEF(3,n): row number of Sn (Scale for comp. n) (thru n+1 output)
      IDEF(4,n): row number of Fn (bioavailability fraction for comp. n)
      IDEF(5,n): row number of Rn (modeled rate for comp. n)
      IDEF(6,n): row number of Dn (modeled duration for comp. n)
      IDEF(7,n): row number of ALAGn (absorption lag for comp. n)

  GG  The array of PK parameters and their eta derivatives.  The  maxi-
      mum  number  of rows in GG is given by IRGG, which is the same as
      constant PG found in file SIZES (See sizes).
      At ICALL = 2:
      GG(k,1,1) contains the value of the kth parameter.
      GG(k,i+1,1) contains its derivative with respect to the ith eta.
      GG(k,i+1,j+1) contains its second derivative with respect to  the
      ith eta and the jth eta (lower-triangular; j=1, ..., i).  (Second
      derivatives are only needed  with  estimation  by  the  Laplacian
      method.)
      At ICALL = 4:
      GG(k,1,1)  contains  the  individual-specific  value  of  the kth
      parameter.  Other columns of GG need not be computed.

 Also see variables in NONMEM modules, NONMEM-PRED modules, and  PREDPP
 modules.
 (See variables in modules)

 REFERENCES: Guide IV, section V.C.5 
 REFERENCES: Guide V, section 7 
 REFERENCES: Guide VI, section III 
