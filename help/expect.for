


 +--------------------------------------------------------------------+
 |                                                                    |
 |                         EXPECTATION BLOCK                          |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Abbreviated code for computation of expectations
 CONTEXT: Abbreviated code

 SAMPLE:
 $ERROR
 IF (ICALL.EQ.5) THEN
  ... expectation block ...
 ENDIF

 DISCUSSION:
 An  expectation block is a block of abbreviated code that is only exe-
 cuted when ICALL=5.  This value of  ICALL  occurs  when  the  marginal
 (MRG_)  data  item is defined in the data set and has a non-zero value
 for some records.  Expectation blocks are not required when  the  mar-
 ginal  data  item is present, but they allow the user additional func-
 tionality.  Such blocks may be present in $PRED, $PK, $ERROR.

 If the MRG_ data item has the value 1 in a data record, PRED- (PK- and
 ERROR-)  defined  items displayed for the record (e.g. in the row of a
 table corresponding to the record) are  expectations.   When  ICALL=5,
 the  expectations are being computed.  With each call to PRED with the
 data record, the value being set in Y is contributing to the  expecta-
 tion  being  computed for the PRED item, and the value being set for a
 PRED-defined item that will be displayed in  a  table  or  scatterplot
 (except  for an item stored in the SAVE region) is contributing to the
 expectation being computed for that item.

 In this example, the displayed value of  PRED  with  a  record  having
 MRG_=1  is  the  expectation of the 0-1 variable taking the value 1 if
 and only if F > 10.

  IF (ICALL.EQ.5) THEN
   Y=0
   IF (F.GT.10) Y=1
   RETURN
  ENDIF

 Special rules apply to expectation blocks and passes through the  data
 with ICALL=5.

 1)   No eta derivatives are computed in an expectation block.

 2)   Calls  to  certain  NONMEM  utility  routines are permitted in an
      expectation block:

      CALL RANDOM(n,R)

      where  n is an integer  1-10.   If  CALL  RANDOM  is  present,  R
      becomes  a  reserved  variable,  and  may not be used outside the
      expectation block.  Multiple calls to RANDOM may be present.

 3)   A RETURN statement may be used in an expectation  block.   If  Y=
      appears  (i.e.  Y  is  assigned  a value) in an expectation block
      before the RETURN (not necessarily the in same  block  containing
      the RETURN), then F is set to Y (F=Y); otherwise F is set to to 0
      (F=0).  If there is no RETURN statement in the expectation block,
      then  as  usual,  F is set to the value of Y assigned by the time
      PRED (PK, ERROR) exits.

 4)   Loops are permitted in an expectation block.  The  syntax  is  as
      follows.

      DO WHILE (condition)
       .. statements ..
      END DO

 5)   If  a mixture model is used, then during calls with ICALL=5, both
      MIXNUM and MIXEST are the index of the subpopulation  into  which
      the individual (whose data record is being passed) has been clas-
      sified.

 (See mrg).

 REFERENCES: None.
