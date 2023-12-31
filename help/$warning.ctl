


 +--------------------------------------------------------------------+
 |                                                                    |
 |                             $WARNINGS                              |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING:  Control  Display  of  NM-TRAN Warning, Data Warning and Data
 Error messages
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $WARNINGS  [NONE] [n] [RESET | NORESET ]
            [WARNINGMAXIMUM= [ NONE | n | (list)] ]
            [DATAMAXIMUM= [ NONE | n | (list)] ]
            [ERRORMAXIMUM=n ]

 SAMPLE:
 $WARNING     WMAX=1,EMAX=9999

 DISCUSSION:
 Optional. Controls the display of NM-TRAN warning  messages  and  data
 error  messages  in  the  current  and  all  subsequent problems until
 another $WARNING record appears.  May appear anywhere in  the  control
 stream  following  $SIZES,  $SUPER (if present) and  $PROBLEM records, |
 but cannot affect the display of warning messages associated with  the
 records that precede it.

 OPTIONS:

 NONE Suppress all warning messages.

 n    Sets  WARNINGMAXIMUM=n, ERRORMAXIMUM=n, DATAMAXIMUM=n.  May be 0,
      in which case it is the same as NONE.

 LIST Print a list of all possible warning messages.  When this  option
      is  present,  the  control stream must consist only of the single
      record:
       $WARNING LIST

 WARNINGMAXIMUM
      Refers to general (non-data) warning messages.  May also be coded
      WARNMAXIMUM  and  shorter  substrings,  e.g.,  WMAX,   WARN.  May
      appear more than once, with cumulative effect.

      WARNINGMAXIMUM=NONE
           Suppress warning messages.

      WARNINGMAXIMUM=n
           Each type of warning message will be displayed no more  than
           n times  (n >= 0).  The default is 20.

      WARNINGMAXIMUM=(list)
           The  list is a series of integers n1 , n2 , n3  ...  ( ni >=
           0 ).  ni gives the the maximum number for the  ith  type  of
           warning message.

           NONE  may  coded for any  ni  if it is the last in the list.
           Sets the maximum number of messages  for  the  ith  and  all
           higher-numbered messages to zero.

      Example:
      $WARNING     WMAX=1
      Sets the maximum for all non-data messages to 1.

 DATAMAXIMUM
      Refers  to data warning messages.  May also be coded DMAXIMUM and
      shorter substrings, e.g. DMAX, DATA.  May appear more than  once,
      with cumulative effect.  Forms and default are identical to those
      for the WARNINGMAXIMUM option.

 ERRORMAXIMUM=n
      Refers to data error messages.  It can also be  coded  ERRMAXIMUM
      or  EMAXIMUM, and shorter strings, e.g., EMAX.  Each type of data
      error message will be displayed no more than n times (n>=0).  The
      default  is  20.   Data  error messages can never be totally sup-
      pressed.  With n=0 and n=1, NM-TRAN displays the first data error
      message and then terminates.

      Example:
       $WARNING EMAX=99999
      Sets  the maximum number of data error messages displayed from 20
      to a much higher value, so that NM-TRAN will not terminate  until
      all data error conditions have been identified.

 RESET
      All counts reset to 0 with each new problem. This is the default.

 NORESET
      Counts  of  how many times a given message has been displayed are
      accumulated over problems.  (A new $WARNING record always  resets
      all counters, but they will not be automatically reset again with
      subsequent problems.)

 A $WARNING record with no options re-establishes  all  default  condi-
 tions as if no $WARNING records had been present.

 REFERENCES: None.
