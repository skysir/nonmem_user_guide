


 +--------------------------------------------------------------------+
 |                                                                    |
 |                   PRED_IGNORE_DATA BLOCK (NM75)                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for ignoring (dropping) data records
 CONTEXT: $PRED, $PK, $INFN abbreviated code

 SAMPLE:
 $INFN
     IF(PRED_IGNORE_DATA_TEST==1) THEN
     PRED_IGNORE_DATA=0
     IF(AGE>35.0) PRED_IGNORE_DATA=1
     IF( ID>10.AND.ID<18.OR.ID>60.AND.ID<70 ) PRED_IGNORE_DATA=1
     RETURN ;Assures no additional computation code in INFN is executed
     ENDIF

 USAGE:
       SUBROUTINE INFN (ICALL,THETA,DATREC,INDXS,NEWIND)
       USE NMPRD_INT, ONLY: PRED_IGNORE_DATA,PRED_IGNORE_DATA_TEST

 DISCUSSION:

 The IGNORE=(list) and ACCEPT=(list) options of $DATA provide a limited
 means of filtering the input data set, which is performed  by  NMTRAN.
 To  provide  more  elaborate  filtering  for  excluding data, PRED can
 request that NONMEM filter out additional data records at  the  begin-
 ning of the run.

 This  is  done  by setting the reserved variable PRED_IGNORE_DATA to a
 non-zero value within $INFN, $PK, or $PRED,  for  each  record  to  be
 ignored.

 It may be useful to package PRED_IGNORE_DATA statements within
  IF(PRED_IGNORE_DATA_TEST==1) THEN
  ...
  RETURN
  ENDIF
 structures to avoid unnecessary code execution.

 If  the  PRED_IGNORE_DATA_TEST or PRED_IGNORE_DATA variables appear in
 abbreviated code, or the option
 $DATA ... PRED_IGNORE_DATA
 is used in the NM-TRAN control stream, then  a  PRED_IGNORE_DATA  pass
 through the NONMEN data file with PRED_IGNORE_DATA_TEST=1 and ICALL=-1
 occurs.  Otherwise it does not.  Therefore, existing code such as  "IF
 (ICALL<=1) THEN ...ENDIF" does not need to be changed.

 The following variables have properly defined values:

 ICALL
 Data record items in DATREC
 NEWIND,NEWL2
 NPROB,IPROB, S1NUM, S2NUM,
 S1NIT,S2NIT, S1IT, S2IT

 No  other variables are properly defined when PRED_IGNORE_DATA_TEST=1.
 For example, the following should not be used:
 THETA, OMEGA, SIGMA, NREP, IREP
 ETA may be used but will be 0 there  are  no  random  variables  in  a
 pred_ignore_data block.

 Typically the NONMEM file that is  input to the pred_ignore_data block
 is FDATA.  FDATA is unaffected by the  pred_ignore_data  block.   How-
 ever,  with  NONMEM  7.5  there  is a new file, FDATA.csv, and records
 excluded by PRED_IGNORE_DATA will not be present in FDATA.csv.

 In a  pred_ignore_data block, the data record has been read by  NONMEM
 and  all  data  record  items  have numeric values.  The (non-Fortran)
 operators .EQN. and .NEN. that can be used with the $DATA  IGNORE  and
 ACCEPT options are not needed and cannot be used in a pred_ignore_data
 block.

 Any other functions of $INFN, such as DATA  item  modification  (i.e.,
 transgeneration  of  the data), RANDOM calls, etc. should be made with
 ICALL==1 or ICALL==0 IF blocks, as before.

 It is possible to restrict PRED_IGNORE_DATA actions  to  a  particular
 problem number:

 IF(IPROB==2.AND.PRED_IGNORE_DATA_TEST==1) THEN
 PRED_IGNORE_DATA=0
 IF(AGE>35.0) PRED_IGNORE_DATA=1
 IF( ID>10.AND.ID<18.OR.ID>60.AND.ID<70 ) PRED_IGNORE_DATA=1
 RETURN
 ENDIF

 RETURN statements may be used.

 EXIT statements may be used.  They act like a RETURN but are otherwise
 ignored.

 REFERENCES: Guide Introduction_7
