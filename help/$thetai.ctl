


 +--------------------------------------------------------------------+
 |                                                                    |
 |                              $THETAI                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives Instructions for Transforming Initial Thetas
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $THETAI  Fortran statements

 SAMPLE:
 $THI
 THETA(1:NTHETA)=LOG(THETAI(1:NTHETA))
 THETAP(1:NTHP)=LOG(THETAPI(1:NTHP))

 THETA=LOG(THETA)

 DISCUSSION:

 The  purpose  of  $THETAI  is  to  transform the initial values in the
 $THETA and $THETAP records.  The record name  may  also  be  coded  as
 $THI.

 In  the  above sample code, it is desired that the thetas by estimated
 within NONMEM in the log domain, but the user wants the convenience of
 inputting and outputting them in the natural domain, such as when per-
 forming linear MU referencing.

 The $THETAI record will convert any initial thetas in a $THETA record,
 or thetas obtained from a chain file, but will not convert thetas from
 an MSF file.  The variance to the theta priors will  be  appropriately
 converted when using $PRIOR NWPRI.

 The  assignment statements may be any Fortran 95 statements.  They are
 copied unchanged to subroutine SUBROUTINE  THETAISUB  in  FSUBS  (also
 found in thetair.f90).

 They  may  include  array  assigment  statements  specifying the whole
 arrays or sections of arrays.

 If the initial estimate for an element of theta is transformed, so  is |
 the upper and lower bounds for that theta, if any.

 Arguments of the subroutine are as follows.

 THETAI
      Values of THETA specified on $THETA records.  Input.

 THETAPI
      Values of THETA specified on $THETAP records (or, if the informa-
      tive names are not used, thetas corresponding to priors, if any).
      Input.

 THETA
      New values of THETA.  Output.

 THETAP
      New values of THETA's for priors.  Output.

 Other reserved variables that may be used are as follows:

 NTHETA
      Number of thetas to be estimated.

 NTHP Number of theta priors.

 NPROB IPROB
      These  can  be  tested  in  IF  statements  so that values may be
      assigned diffently for different problems.

 If the range is not specified, NONMEM to supply the range (which is by
 default NTHETA+NTHP).

 When appropriate, for reporting thetas, the inverse function should be
 supplied, e.g., with the samples above:

 $THR
 THETAR(1:NTHETA)=EXP(THETA(1:NTHETA))
 THETAPR(1:NTHP)=EXP(THETAP(1:NTHP))
 or
 THETAR=EXP(THETA)

 Note that the assignment occurs after the NONMEM  control  stream  has
 been processed, so that errors in assignment of THETA's are not found.
 E.g.
 $THI
 THETA=0.
 This will set all theta's to 0, and there will be  no  specific  error
 message from NONMEM, though most likely the run will fail.

 If  initial estimates of all or part of THETA are omitted, NONMEM per-
 forms a search for  missing  initial  estimates  as  its  first  task.
 (These  are  printed  on  a  separate  page under the heading "INITIAL
 PARAMETER ESTIMATE").  This search occurs before the transformation by
 $THETAI.

 Values  of  INITIAL  ESTIMATE OF THETA in NONMEM output are those from
 the $THETA/$THETAP records.  Values of THETA in the root.ext files are
 as set in $THI.

 $THI may be used with $THR, but not necessarily.

 Another example is rescaling thetas.  E.g., suppose in CONTROL5

 $THETA  (.1,3,5) (.008,.08,.5) (.004,.04,.9)

 is replaced with

 $THETAI
 THETA=THETAI/10.
 $THETA  (1,30,50) (.08,.8,5) (.04,.4,9)

 The  results will be identical, except for the values of INITIAL ESTI-
 MATE OF THETA in the NONMEM report. The values in  .ext  will  not  be
 affected.

 REFERENCES: Guide Introduction_7
