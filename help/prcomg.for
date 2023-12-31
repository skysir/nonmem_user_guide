


 +--------------------------------------------------------------------+
 |                                                                    |
 |                      DIFF EQ SOLVER SETTINGS                       |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: PREDPP global variables
 CONTEXT: For use with PREDPP

 This feature is not fully documented.  The interested user may be able
 to obtain more information by studying  the  appropriate  sections  of
 PREDPP code.

 USAGE:
 $PK
 "FIRST
 " USE PRCOM_INT, ONLY: METH,MITER,IMAX,ISTFLG,INTFLG
 "MAIN
 " IMAX=200000

 DISCUSSION:
 These  variables  allow the user to over-ride certain default settings
 in  ADVAN6,  ADVAN8,  ADVAN9,  ADVAN13,  ADVAN14,  ADVAN15,   ADVAN16,
 ADVAN17, ADVAN18, SS6, and SS9.

 DES_DER, MITER, and METH

      PREDPP sets MITER and METH to default values with every new indi-
      vidual or reset record.  METH is  the  solver  method  type,  and
      MITER determines whether analytical Jacobian
       is  to  be  used (MITER=1, based on DA() array calculated in DES
      subroutine), or Jacobian is to be numerically determined  by  the
      solver  (MITER=2).  As  of  NM75,  MITER  may  be accessed by the
      reserved variable DES_DER,  without  requiring  the  "FIRST"  and
      "USE" header lines, or verbatim code. For example:
           $PK
           DES_DER=1

      ADVAN8
           METH=2 is always used.
           If user set MITER, use his value.  Otherwise, use MITER=2.

      ADVAN9
           If  METH=1, then Adams implicit method, and METH=2 (default,
           BDF method).

           If user did set MITER, use his value.

      ADVAN13
           METH set in $PK is not used. LSODA determines whether  ADAMS
           NON-stiff (METH=1) or BDF stiff (METH=2) is to be used as it
           works its way through the problem.
           MITER (or DES_DER) should be set only to 1 (analytical Jaco-
           bian) or 2 (numerical Jacobian determined by LSODA). Usually
           should be left alone.

      ADVAN14
           METH=2 by default. Could be set in cvodeu.f90 or $PK.
           Setting MITER=1 (or DES_DER=1) is required if it is  desired
           that  the  Jacobian be analytically evaluated for ODE models
           (ADVAN8, 9, 13, 14, 15, 16, 17), for IMP, SAEM, BAYES  prob-
           lems  that normally do not have first derivatives turned on.
           In addition to setting DES_DER=1 in $PK, various options  of
           Jacobians  (analytical/numerical,  full/band, etc) should be
           set in CVODEU.f90, where there are  many  other  things  the
           user can set.

      ADVAN15
           METH is not used. It is always METH=2 (BDF stiff).
           Setting  MITER=1 is required if it is desired that the Jaco-
           bian be analytically evaluated for ODE  models  (ADVAN8,  9,
           13,  14, 15), for IMP, SAEM, BAYES problems that normally do
           not have first derivatives turned on.  In addition  to  set-
           ting DES_DER=1 in $PK, various options of Jacobians (analyt-
           ical/numerical, full/band, etc) should be set  in  IDAU.f90,
           where there are many other things the user can set.

      ADVAN16
           METH  and  MITER  are  not used.  Instead, IJAC would set in
           RADAR5u.f90, but presently this is not functional

      ADVAN17
           METH and MITER are not used.  Instead,  IJAC  would  set  in
           RADAR5u.f90, but presently this is not functional

      ADVAN18
           METH and MITER are not used.

 IMAX

      ADVAN6, ADVAN8, ADVAN9, ADVAN13, ADVAN14, ADVAN15, ADVAN16,
           ADVAN17, ADVAN18, SS6
           The variable MAXCAL gives the maximum  number  of  calls  to
           FCN1  (ADVAN6,  ADVAN8,  ADVAN13, ADVAN14, ADVAN16, ADVAN18,
           SS6) or RES (ADVAN9, ADVAN15, ADVAN17) during an integration
           interval.
           Each of the above routines sets MAXCAL to the value given by
           MAXFCN (a parameter in the MODULE SIZES)  at  the  start  of
           each  integration  interval unless the user supplies a value
           in IMAX, in which case the user's value is used.

 ISTFLG

      ADVAN9 ADVAN15, ADVAN17
           ISTFLG controls how ADVAN9 calls LSODI1 and what it does  if
           LSODI1  returns  and  indicates  that an integration failed.
           ISTFLG is set to 0 (default) at ICALL=0.  If changed by  the
           user, it retains the changed value until the user changes it
           again.  ISTATE is a variable that is passed from  the  ADVAN
           to  LSODI1.   ISTATE=1  indicates  that  the  integration is
           starting. ISTATE=2 indicates that this call is  a  continua-
           tion  from  a  prior  successful  integration.  Default: Use
           ISTATE=1 for the first integration, and ISTATE=2 for a  con-
           tinuation  when  nothing  external to the ADVAN has changed.
           In case of failure with ISTATE=2,  restore  original  inputs
           and try again with ISTATE=1.

      ISTFLG=1
           Never try ISTATE=2, always use ISTATE=1.

      ISTFLG=2
           Never retry (only try ISTATE=2).

      ADVAN15 also uses ISTFLG for calls to IDA, but uncertainty exists
      as to what its usefulness is.

 INTFLG

      ADVAN6, SS6, ADVAN8, ADVAN13, ADVAN14, ADVAN16, ADVAN18

      INTFLG stands for "Integration Flag" and affects  the  number  of
      calls  to  the integrating subroutine during each advance.  It is
      only of interest when second derivatives of the state vector  are
      calculated.

      It  is  present  because  there may be some trade-off between run
      time and accuracy of computation.  More calls to the  integrator,
      with  a  smaller  number  of derivatives obtained with each call,
      result in longer run times, but might also produce more  accurate
      derivatives.  Also, it might provide more consistent computations
      when the number of compartments and/or etas is to be changed.

      Default is -1.  User may set to any other value  in  user-written
      code (e.g., with verbatim code in $PK).

      PREDPP  examines  this  value  after  the first call to PK for an
      individual record, so that it can be  set  on  an  individual-by-
      individual  basis.   PREDPP  examines  INTFLG  when NEWIND=0 or 1
      (i.e., at the start of an individual's data).   Presumably,  this
      is  the  only  time  that  NONMEM  might  change  the LVOUT array
      (See non-active eta list for pred).

      When  INTFLG  is  set  to  any  other  value   than   -1,   ADVAN
      6,8,13,14,16,18   and SS6 calculate second derivatives "one group
      at a time".

      E.g., with ADVAN6, suppose etas 1, 2, 3 are active.

           1) call DVERK to obtain 2nd derivatives 1,1
           2) call DVERK to obtain 2nd derivatives 1,2 and 2,2
           3) call DVERK to obtain 2nd derivatives 1,3 and 2,3 and 3,3

           (Each calculation involves the integration of the state vec-
           tor,  augmented  by the relevant first derivative(s) and the
           second derivatives for one eta.)

           Thus, the maximum  number  of  differential  equations  that
           would ever be integrated at one time is

               PW=2*PE*PM+PM
           Where: (PE is the maximum number of etas; currently 10)
           (PM is the maximum number of compartments - 1; currently 9)
           PW is the size defined for various work arrays in the source
           code.
           It is currently defined as above (189).

      When INTFLG is -1, the ADVAN routines  and  SS6  makes  the  most
      efficient  use  of  the work arrays when computing second deriva-
      tives, to reduce the number of calls to DVERK.

      At NEWIND=0 or 1, PREDPP looks at the number of active  etas  for
      this  individual  and  the  number  of  user-defined compartments
      defined by the MODEL subroutine (NCM), and creates  a  scheme  to
      calculate  as  many groups of 2nd. derivatives at once as it can.
      (Although compartments may be turned on and off within  the  data
      set, PREDPP does not revise the scheme of integration each time.)

           E.g.,  with  the  current values of PE, PM, and PW, here are
           some schemes of integration.  Under "nth. call" are the etas
           whose  second  derivatives  are  computed  with that call to
           DVERK.

           # compts.  1st. call     2nd. call    3d. call  4th. call
              9       etas 1 - 5    etas 6 - 7   eta 8     eta 9
              8       etas 1 - 5    etas 6 - 7   eta 8     eta 9
              7       etas 1 - 5    etas 6 - 7   eta 8 - 9
              6       etas 1 - 6    etas 7 - 8   eta 9
              5       etas 1 - 7    etas 8 - 9
              4       etas 1 - 8    etas 9
              3       etas 1 - 9

      Changing the size of the work arrays:

      If the system is large enough that integration involves more than
      one call to DVERK, and the user would like all second derivatives
      to be computed with a single call to DVERK, the source code  must
      be changed to define a larger value for PW.

      To  integrate  the maximum number of etas and compartments in one
      call, set:

        PW=PM*(1+PE+PE*(PE+1)/2)

      With the current values of PE and PM, PW=594

      Suppose PW is changed, but by accident is made smaller  than  the
      default  (2*PE*PM+PM=189).   Then for problems with large numbers
      of compartments and/or etas, the work arrays will  not  be  large
      enough, with either INTFLG=-1 or INTFLG!=-1.

      A new error message exists in PREDPP, for which PRED exit code is
      2 (always abort).

      WORK ARRAYS ARE TOO SMALL  FOR  2ND.  DERIVS.   INCREASE  PW,  OR
      DECREASE NO. OF. COMPTS AND/OR ETAS, OR USE DERIV2=NO

      Again, this message cannot occur unless the source code of PREDPP
      is changed incorrectly.

 Location prior to NONMEM 7: prcomg

 REFERENCES: None.
