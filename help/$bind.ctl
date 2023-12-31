


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $BIND                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Define data values used by $PK, $DES, and $AES
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $BIND  [value1]  [value2] ...

 DISCUSSION:
 $BIND  is  optional.   It  may  be  used  when $PK abbreviated code is
 present and this code requests that the PK  subroutine  be  called  at
 additional  or  lagged dose times (CALLFL=-2).  It is used to override
 the default values of user-defined variables used in the code when the
 PK  routine is called at these particular dose times.  It is also used
 to override the default values of user-defined variables used in  $DES
 and $AES abbreviated code during an advance to an additional or lagged
 dose time.

 $BIND has no effect when PK is called at a model event time (MTIME).   |

 Let t be a time at which an additional or lagged dose enters the  sys-
 tem.   If t1 is the greatest event time not exceeding the time t, then
 the "last event record" is the last event record with event  time  t1,
 and  the  "next  event  record" is the event record following the last
 event record.
 (The term "last" is similar to the word "previous" in this context.)
 The event time on the next event record will exceed time t.

 OPTIONS:

 The positions of the values correspond to the positions of data  items
 in the $INPUT record in a 1-to-1 manner.  Each value is one of:
   DOSE (Use the value from the dose record.)
   NEXT (Use the value from the next event record.)
   LAST (Use the value from the last event record.)
   SKIP (Ignore this data item.)
   DROP (Ignore this data item.)
   -   (Default.)

 For  user  (concomitant)  data  items,  the default is NEXT but any of
 DOSE, NEXT, LAST may be specified.

 For the PREDPP data item TIME, the default is NEXT.  Only  -  or  NEXT
 may be specified.

 For  all other PREDPP or NONMEM data items, the default is DOSE.  Only
  -  or DOSE may be specified.

 Option DROP in the $BIND record is optional.  It  is  ignored  by  NM-
 TRAN,  but helps the user maintain a 1-to-1 relationship between posi-
 tions in $BIND and positions in $INPUT.

 A $BIND record with all defaults specified, such as
 $BIND - - - - - - -
 has the same effect as if no $BIND record were present.

 EXAMPLE:

   $INPUT ID TIME DATE=DROP AMT DV WGT  PREP  X   HGT
   $BIND   -   -    DROP     -   - NEXT DOSE LAST  -

     ID, AMT, DV  have the values from the initiating dose record.
     TIME has the value from the next event record.
     WGT and HGT have values from the next event record.
     PREP has the value from the initiating dose record.
     X has the value from the last event record.  This record  will  be
     the dose record if there is no other event record between the dose
     record and the next event record.

 $INPUT and $BIND records can be interleaved to help maintain a  visual
 relationship. The above example could have been coded:

   $INPUT ID TIME DATE=DROP
   $BIND   -   -    DROP
   $INPUT AMT DV  WGT PREP X    HGT
   $BIND  -   -  NEXT DOSE LAST  -

 $BIND  may not specify a position beyond the last position defined via
 $INPUT.  It may specify fewer positions, in which case defaults  apply
 to the remaining data items.

 Changes  to  $BIND, like changes to $INPUT, cause changes to generated
 code.  Thus, an existing NONMEM executable cannot be re-used when  the
 $BIND and/or $INPUT records are changed.

 The $BIND record only applies under the following circumstances:

 (1)  There exists a dose that is either one or both of the two follow-
      ing types of doses:

      a)   An  additional  dose,  subsequent  to  the  initiating  dose
           (ADDL>0, II>0 in the dose event record).  Such a dose enters
           the system at the "additional dose time."

      b)   A lagged dose.  (with the corresponding dose  event  record,
           PK  specifies ALAGi>0, where i is the index of the dose com-
           partment).  Such a dose enters the  system  at  the  "lagged
           dose time."

 (2)  The  PK  subroutine  computes parameters that depend on values in
      the data record which are not constant for the individual,  i.e.,
      parameters depend on time-varying data items.

 (3)  The PK subroutine is called when the dose enters the system.
      That  is, $PK contains the pseudo-statement CALLFL=-2, requesting
      that the PK routine be called to compute values of the PK parame-
      ters at additional and lagged dose times.  When $PK does not con-
      tain this pseudo-statement, there is no  such  call,  and  at  an
      additional or lagged dose time, the PK parameters have those val-
      ues computed with the next event record.

 When the $PK code is implemented with an event  record,  the  variable
 DOSTIM  is  0.  When it is implemented at an additional or lagged dose
 time, the value of this variable is the time in question.  By default,
 data items used in abbreviated code have values from either the initi-
 ating dose record (DOSREC) or the next event record (EVTREC),  accord-
 ing to this rule:

 NONMEM  data items and PREDPP data items (other than TIME) have values
 from the original dose event record (DOSREC).  These are:

 DV MDV ID L2 MRG_ RAW_ REPL_ EVID AMT RATE SS II CMT  PCMT  CALL  CONT
 ADDL DATE DAT1 DAT2 DAT3

 TIME  and  all  user-defined  data items have the values from the next
 event record (EVTREC).  This default can be overridden using the $BIND
 record.

 EXAMPLE:

 The  $BIND  record  is  a  convenience; it does nothing that cannot be
 accomplished with abbreviated code.  Suppose X is a user-defined  data
 item.  The following three fragments of code create a variable XB with
 the same values at calls with additional and lagged doses that X would
 have if there were a $BIND record specifying the following for X.

 DOSE:
   IF (DOSTIM.EQ.0) THEN
   XB=X
   ELSE
   XB=DOSREC(X)
   ENDIF

 NEXT:
   XB=X

 LAST:
   IF (DOSTIM.EQ.0) THEN
   XB=X
   ENDIF

   In this case, when DOSTIM>0, XB retains its value from the previous
   call to PK.

 (See bind example).

 REFERENCES: Guide IV, section V.C.2 , V.C.5 
