


 +--------------------------------------------------------------------+
 |                                                                    |
 |                        INDICATOR VARIABLES                         |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Coding technique
 CONTEXT: Abbreviated code

 DISCUSSION:
 An indicator variable is a variable whose value is 0 or 1.   It may be
 identified with an input data item, or it may be a variable defined in
 the abbreviated code.  It is used to make a choice in a computation.

 Suppose,  for  example, ICU is a variable which is either 0 or 1.  The
 following code can be used:

 IF (ICU.EQ.0) THEN
 CLM=TVCLM+ETA(1)
 ELSE
 CLM=TVCLM+ETA(2)
 ENDIF

 This can be coded unconditionally using an indicator variable:
 CLM=TVCLM+(1-ICU)*ETA(1)+ICU*ETA(2)

 Unconditional code is preferred when MU variables are computed.        |

 In Guide V, Chapter 12, appears this example involving observations of
 two  types,  CP  and effect.  The latter is modeled by the Emax model.
 Observations of both types are recorded in the DV data item.   A  data
 item  called  TYPE  with  values  1  and  2  is  used  to  distinguish
 between  them.  This data item may be used to compute the  appropriate
 prediction in the $ERROR block, as follows:

 $ERROR
  EMAX=THETA(5)+ETA(3)
  C50=THETA(6)+ETA(4)
  EFF=EMAX*F/(C50+F)
  Y1=EFF+ERR(1)
  Y2=F+ERR(2)
  IF (TYPE.EQ.2) THEN
  Y=Y2
  ELSE
  Y=Y1
  ENDIF

 The last five lines of this example can also be coded as follows:

  Q=1
  IF (TYPE.EQ.2) Q=0
  Y=Q*Y1+(1-Q)*Y2

 A  more  general  technique  is  to use two indicator variables.  This
 technique can easily be extended to the case of three or more choices.

 $ERROR
  EMAX=THETA(5)+ETA(3)
  C50=THETA(6)+ETA(4)
  EFF=EMAX*F/(C50+F)
  Y1=EFF+ERR(1)
  Y2=F+ERR(2)
  Q1=0
  Q2=0
  IF (TYPE.EQ.1) Q1=1
  IF (TYPE.EQ.2) Q2=1
  Y=Q1*Y1+Q2*Y2

 REFERENCES: Guide V, section 7.5.3 , 7.5.4 ,  12.5  
