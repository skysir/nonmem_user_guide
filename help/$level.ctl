


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $LEVEL                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Specifies Nested Random Levels Above Subject ID               |
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $LEVEL item=(n1[m1] , n2[m2] ... ) ...

 SAMPLE:
 $LEVEL

 DISCUSSION:

 Identifies one or more Super ID data items.

 item Item  is the name of a data item listed on $INPUT.  It defines an
      additional nested random level.  and is referred to as  a  "super
      ID"  data  item.   The  first super ID item defines an additional
      random nesting level above that of subject  ID.   More  than  one
      super  ID item may be listed on $LEVEL.  Each subsequent super ID
      item defines an additional nesting level above that of the previ-
      ous nesting level, i.e., above the previous super ID.

 nk[mk]
      States  that  ETA(nk)  is associated with this super ID item, and
      ETA(mk) is nested within ETA(nk).

 With NONMEM 7.4, a short-hand notation  may  be  used  to  describe  a |
 series  of values of nk.  A sequence of values for nk can be described |
 as                                                                     |
  start TO end BY interval                                              |

      TO is required.  The character : may be used instead of  TO.   BY |
      is optional. Default is 1.  The value of BY may be negative.      |

 If  the  second  value of BY differs, the same syntax may also be used |
 for mk.                                                                |

 Nesting levels below the subject ID is modelled as with previous  ver-
 sions of NONMEM.
 (See Interoccasion_variability example).

 The  order  that super ID's are listed on $LEVEL defines their nesting
 level.  The order that standard and super ID's are  listed  on  $INPUT
 (i.e.,  the order in which they appear in each record of the data set)
 is immaterial.

 When $LEVEL is used with FOCE ($ESTM METHOD=1),  the  SLOW  option  is
 required, and MATRIX=R is required with $COV.

 EXAMPLE:

 $INPUT ... ID ... SID ... CID ...
  ...
 $PK
 MU_1=THETA(1)
 MU_2=THETA(2)
 CL=DEXP(MU_1+ETA(1)+ETA(5)+ETA(9))
 V1=DEXP(MU_2+ETA(2)+ETA(6)+ETA(10))
  ...
 $LEVEL
 SID=(5[1],6[2])
 CID=(9[5],10[6])

 The  data  item  named SID is the site ID.  The data item named CID is
 the country ID.  There are several sites  belonging  to  one  country,
 some  other  sites  belonging to another country, etc.  For clearance,
 eta(9) is the country variability that has nested in it the site vari-
 ability eta(5), which in turn has nested in it the subject variability
 (the standard ID data) eta(1).  For V1, eta(10) is the  country  vari-
 ability  that  has  nested in it the site variability eta(6), which in
 turn has nested in it the subject variability (the standard  ID  data)
 eta(2).

 An alternate way of coding the $LEVEL records is                       |

 SID=(5 to 6[1])                                                        |
 CID=(9 to 10[5])                                                       |

 NONMEM  performs  appropriate summary statistics for eta(5), and makes
 the appropriate constraints on eta(5), so eta(5) changes by site, that
 is, by every SID value change, and not by every ID value change.

 The  above  method, using $LEVEL, is a linearized approximation at the
 super ID level, and takes advantage  of  a  dual  OBJ  function  call,
 freely  allowing all etas to vary on the first call of OBJ, then aver-
 aging the SID etas, fixing them to these averages, and  going  through
 another  OBJ call to allow the subject (ID) etas to be assessed.  This
 approximation method works very well for the EM and Monte Carlo  meth-
 ods, and reasonably well for the FOCE/Laplace methods.

 To perform an exact analysis, separate thetas must be defined for each
 value pertaining to a super ID data item, so that theta is shared only
 by the subjects with the particular SID value.  $LEVEL is not used.
 (See superid3_6).

 If  there  are  multiple  $PROBLEM records, $LEVEL  should be restated |
 with each problem for which it is still relevant.  For  example,  this |
 is the case with $SUPER problems and $LEVEL.

 See  also  LEVWT option of the $ESTIMATION record (NM74).  By default,
 LEVWT=0, and weights each level value equally, regardless of number of
 subjects per level value. If you wish to weight according to number of
 subjects for that value, set LEVWT=1 on the $EST record.

 REFERENCES: Guide Introduction_7
