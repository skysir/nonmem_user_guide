


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      COM COMACT COMSAV COMRES                      |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Variables related to the module NMPRD4
 CONTEXT:  Abbreviated code, verbatim code, user-supplied routines, NM-
 TRAN

 Reserved labels starting COM are related  to  each  other.   They  all
 refer  in  some  manner to the module NMPRD4.  (This module was origi- |
 nally a global COMMON; hence the use of the letters "COM".)

 COM(n) (in $TABLE and $SCATTER records)
      Certain positions of MODULE NMPRD4 may be reserved, and the vari-
      ables  stored in these positions may be displayed by listing them
      as COM(1), COM(2), etc.  on $TABLE and $SCATTER records, e.g.,
       $TABLE COM(1) COM(2)

 COM(n) (in abbreviated and verbatim code)
      Certain positions of MODULE NMPRD4 may be reserved, and the vari-
      ables stored in these positions are referenced in abbreviated and
      verbatim code as COM(1), COM(2), etc.

 COMACT (in any user code)
      Reserved variable COMACT is set by NONMEM. It may  be  tested  in
      PRED  (e.g., in abbreviated code, verbatim code, or in user-writ-
      ten routines) to determine when NONMEM is making a copying  pass,
      i.e., when the data records are being passed to PRED for the pur-
      pose of computing values of  variables  which  will  be  obtained
      (i.e.  copied)  from  NMPRD4 for tables and scatterplots.  NONMEM
      only makes a copying pass when PRED-defined items are  listed  in
      $TABLE  or  $SCATTER records.  There may be a few copying passes.
      With the first copying pass, the value of COMACT is 1.  As  copy-
      ing  passes  proceed,  the value of COMACT may remain the same or
      increase.  The values used in tables and scatterplots  are  those
      copied from NMPRD4 with the last copying pass.
      COMACT=0: This is not a copying pass.
      COMACT=1:  This is a copying pass with final thetas and zero-val-
      ued etas.
      COMACT=2: This is a copying pass with  final  thetas  and  condi-
      tional estimates of etas.
      COMACT=3: This is a copying pass with conditional (nonparametric)
      estimates of etas.
      Such a pass takes place when the  control  stream  includes  this
      record:
       $NONPARAMETRIC ETAS

 COMRES=n1 (option of $ABBREVIATED record)
      COMRES  ("common  reserve")  is  an  option  of  the $ABBREVIATED
      record.  It gives instructions to NM-TRAN about NMPRD4.
      COMRES=-1:  Do not store any variables in module NMPRD4.
      COMRES=0:  Store variables in NMPRD4,  but  do  not  reserve  any
      positions (the default).
      COMRES=n1:   Store  variables in NMPRD4, and reserve the first n1
      positions.

 COMRES=-1 (in abbreviated code)
      The pseudo-assignment statement COMRES=-1  may  be  used  in  any
      block  of abbreviated code to prevent any variable defined in the
      block from being stored into NMPRD4.

 COMSAV=n2 (option of $ABBREVIATED record)
      Values of variables displayed  in  tables  and  scatterplots  are
      obtained  from  module  NMPRD4.   There are particular times when
      data records are passed to PRED  for  the  purpose  of  obtaining
      these  values;  these are called copying passes.  The SAVE region
      of module NMPRD4 is the initial part of NMPRD4.  If a variable is
      stored  in  the  SAVE region, then the value of the variable com-
      puted with a given data record during  a  copying  pass  will  be
      found  in  NMPRD4  when the same record is passed during the next
      copying pass, i.e. it will have  been  saved  from  the  previous
      copying  pass.  This is in contrast to the usual behaviour, where
      with a given data record, the value in NMPRD4 is the  value  com-
      puted with the previous data record.
      n2  is  the  initial  size of the SAVE region, i.e. the number of
      positions in this region.  n2 =0 is the default  value.   n2  may
      not exceed n1.
      The  SAVE region has size n2 initially, but NM-TRAN may extend it
      if SAVE variables are used.  However, if n2 =-1, the SAVE  region
      is  not  to  be extended, and there is to be no SAVE region alto-
      gether.
      (See copying_block).

      NM-TRAN causes the generated routines to store the value of  COM-
      SAV at ICALL=1.

 COMSAV=n2 (in user written code)
      In  the absence of abbreviated code, COMSAV may be set by a user-
      written PRED routine at ICALL<=1, or at COMACT=1  with  the  very
      first data record.  n2 is as described above.

 (See COMACT,COMSAV)
 (See PRED-Defined Variables)
 (See abbreviated, abbreviated_code, displayed_values).

 REFERENCES: Guide IV, section III.B.7 , IV.E.2 
 REFERENCES: Guide VI, section III.J 
