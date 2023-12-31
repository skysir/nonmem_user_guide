


 +--------------------------------------------------------------------+
 |                                                                    |
 |                   CORRELATION ACROSS L2 RECORDS                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM-PRED global variables
 CONTEXT: PRED and ERROR routines

 USAGE:
      CORRL2(n,m)=...

 DISCUSSION:

 An  individual  record  is  divided into L2 records.  An L2 record may
 contain one or more observations (on one or more separate data records
 respectively),  in  which  case it is called an observation-L2 record.
 The values of epsilons used in the intraindividual model may be corre-
 lated across the observations contained in the L2 record, and thus the
 L2 record may define a multivariate observation - the L2  observation.
 (When  all  L2  observations  in  the data set are univariate, L2 data
 items need not appear, and when L2 data items do  not  appear,  NONMEM
 assumes that each data record is a distinct L2 record.)

 By  default,  the values of a given epsilon are statistically indepen-
 dent across  L2  observations  within  an  individual  record.   Using
 CORRL2,  however, these values may be correlated.  More precisely, the
 values of the epsilons associated with the mth diagonal block of SIGMA
 may  be  correlated  across L2 observations, and it will be understood
 that for two different epsilons (eps1 and eps2, say)  associated  with
 the  mth  block,  the  correlation  between  the values of eps1 for L2
 observations A and B will be taken to be the same as  the  correlation
 between the values of eps2 for these same two L2 observations.

 With NONMEM 7.3, reserved variable CORRL2 is used and code such as the |
 following may be used in $ERROR and $PRED blocks.

 Proceed as follows.  With the first data record of  the  nth  observa-
 tion-L2 record, and with respect to the values of the epsilons associ-
 ated with the mth diagonal block of SIGMA, the PRED routine should set
 CORRL2(k,m), for k=1,...,n, to the correlation  between the values for
 the kth L2 observation  and  the  nth  L2  observation.   (CORRL2(n,m)
 should be set to 1.0; in particular, with the first data record of the
 1st observation-L2 record, CORRL2(1,m) should be set to 1.0.)

 E.g. Suppose that the L2 observations are  bivariate  and  chronologi-
 cally  ordered by a TIME data item, and suppose that the intraindivid-
 ual model has two epsilons (one for each element of the  bivariate  L2
 observation),  each  associated with the same diagonal block of SIGMA,
 Then the values of these epsilons may be autocorrelated across the  L2
 observations,  as specified by the above code. But if the two epsilons
 are associated with two different diagonal blocks, then one might  use
 this  code,  in  which  each  sigma  block has its own autocorrelation
 parameter theta(4) or theta(5):

 $ABBREVIATED DECLARE T1(NO)
 $ABBREVIATED DECLARE INTEGER I1
 $ABBREVIATED DECLARE DOWHILE J1

 IF (NEWL2==1.AND.EVID==0) THEN
   I1=I1+1
   T1(I1)=TIME
   J1=1
   DO WHILE (J1<=I1)
   CORRL2(J1,1)=EXP(-THETA(4)*(TIME-T1(J1))) ; Sigma block 1
   CORRL2(J1,2)=EXP(-THETA(5)*(TIME-T1(J1))) ; Sigma block 2
   J1=J1+1
   ENDDO
 ENDIF

 $THETA 2  ;[CL]
 $THETA 30 ;[V]
 $THETA 0.05 ;[Rho1]
 $THETA 0.075 ;[Rho2]

 If you wish to have more control as to when  the  CORRL2  is  used  to
 model  among  the  individual data records within the L2 records, con-
 sider the following example:

 $ABBR DECLARE T1(NO),
 $ABBR DECLARE INTEGER I1, DOWHILE J1

 $ERROR
 IF (NEWIND.NE.2) I1=0

 IF (NEWL2==1) THEN
   I1=I1+1
   T1(I1)=TIME
   J1=1
   DO WHILE (J1<=I1)
   IF(MV1==0.0) CORRL2(J1,1)=EXP(-THETA(4)*(TIME-T1(J1)))
   IF(MV2==0.0) CORRL2(J1,2)=EXP(-THETA(5)*(TIME-T1(J1)))
   J1=J1+1
   ENDDO
 ENDIF
 IF(CMT==1) Y=F+F*EPS(1)+EPS(2)
 IF(CMT==2) Y=F+F*EPS(3)+EPS(4)
  ...
 $SIGMA BLOCK(2)
 0.3
 0.001 0.04
 $SIGMA BLOCK(2)
 0.7
 0.001 0.08

 Here MV1 and MV2  data item that is assumed to exist in the data  set,
 and MV signals the desire to use CORRL2 on an L2 observation within an
 L2 record,  on the first data record of the L2 record.  The  following
 is a section of data to indicate more clearly how this may be set up:

 ID   TIME L2 MV1 MV2 CMT DV
 1    2.0  1  0   0   1   0.8
 1    2.0  1  0   0   2   10.0
 1    4.0  2  0   0   1   0.6
 1    4.0  2  0   0   2   12.0
 1    6.0  3  1   1   1   0.3
 1    6.0  3  1   1   2   20.0
 1   10.0  4  1   0   1   0.1
 1   10.0  4  1   0   2   15.0
 1   14.0  5  1   1   1   0.03

 Only  data  points  with times 2, 4, and 10 hours will incorporate the
 correlation of the CORRL2 model. Furthermore the  record  of  time  10
 hours has MV1=1 and MV2=0, so the observation record in the 10 hour L2
 record that uses Sigma block 2 will have a CORRL2 assessment,  whereas
 the observtation in the 10 hour L2 record that uses Sigma block 1 will
 not be CORRL2 modeled.  In the above example, Sigma block  1  consists
 of EPS(1) and EPS(2), and Sigma block 2 consists of EPS(3) and EPS(4),
 with the CMT value determining which  data  records  use  which  Sigma
 blocks (epsilons sets) based on how Y is evaluated.

 RESTRICTIONS:

 With  versions  of  NONMEM  before  7.3,  C should be used rather than
 CORRL2.  Because C is not recognized by NM-TRAN, and because of  other
 restrictions  regarding  abbreviated  code,  a  specification of C, as
 above, within a block of abbreviated code, must be done using verbatim
 code.  See help for that version of NONMEM.

 SIMULATION WITH POPULATION DATA AND AUTO-CORRELATION                   |

 If  population  data are simulated, the correlations must be stored in |
 CORRL2 before the NONMEM utility routine SIMEPS is called.  With  NON- |
 MEM  7.3,  code  such as the following may be used in $ERROR and $PRED |
 blocks.  Because the value of CORRL2 is set before SIMEPS  is  called, |
 NM-TRAN omits the usual default call to SIMEPS.

 $ABBR DECLARE T(NO), INTEGER I, DOWHILE J

 IF (ICALL == 4) THEN
 IF (NEWIND /= 2) I=0
 IF (MDV == 0) THEN
    I=I+1
    T(I)=TIME
    J=1
    DO WHILE (J <= I)
    CORRL2(J,1)=EXP(-THETA(4)*(TIME-T(J)))
    J=J+1
    ENDDO
 ENDIF
 CALL SIMEPS(EPS)
 ENDIF

 RESTRICTIONS:

 With  versions  of  NONMEM  before  7.3,  C should be used rather than
 CORRL2.  Because C is not recognized by NM-TRAN, and because of  other
 restrictions  regarding  abbreviated  code,  a  specification of C, as
 above, within a block of abbreviated code, must be done using verbatim
 code.  Since NM-TRAN generated code has a call to SIMEPS in its second
 section (see Guide IV), this means that correlations must be  computed
 and  stored  using  verbatim code in the "FIRST block. (Details can be
 supplied on request.)

 SINGLE-SUBJECT DATA

 "Single-subject data" with correlated residual error, can be simulated
 and  analyzed.   To  do  this, though, a technique is needed which can
 always be used with such data: the data are handled  as  data  from  a
 population  sample  with a single individual, and OMEGA is constrained
 to be 0.

 (See Simulation:_SIMEPS_Error_Code)

 Location prior to NONMEM 7: nmprd5

 REFERENCES: None.
