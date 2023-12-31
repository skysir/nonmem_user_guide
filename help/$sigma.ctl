


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $SIGMA                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Supplies initial estimates for the NONMEM SIGMA Matrix
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SIGMA  [DIAGONAL(n) | BLOCK(n) | BLOCK(n) SAME(m) | BLOCK SAME(m)]
         [[value1]  [value2]  [value3] ...
         [(value,value...) xn]
         [BLOCK(n) VALUES(diag,odiag)]
         [label=value] ...
         [BLOCK(n) [NAMES (label1,...,labeln)] [VALUES (diag,odiag)]
         [FIXED] [UNINT]
         [VARIANCE|STANDARD] [COVARIANCE|CORRELATON] [CHOLESKY]

 SAMPLE:
 $SIGMA BLOCK(3)  6. .005 .3 .0002 .006 .4

 DISCUSSION:

 Gives initial estimates and constraints for elements of one or several
 blocks of the SIGMA matrix, i.e., the matrix of variances and  covari-
 ances  of  the  eps  variables  in the statistical model.  This record
 should appear only if the statistical model  contains  eps  variables.
 Multiple  $SIGMA  records  may  be  used  to define multiple blocks of
 SIGMA.  The order of the appearance of all  blocks  over  all  records
 corresponds to the order of the blocks in SIGMA.

 If the initial estimates are omitted for any element(s) of SIGMA, then
 NONMEM will try to obtain the initial estimates.

 OPTIONS:

 There are six forms:

 1.   $SIGMA   [DIAGONAL(n)] [ v11 v22 v33 ... vnn ]

      This gives the initial estimates of the diagonal  elements  of  a
      diagonal block of SIGMA.  E.g.,
        $SIGMA .04 .12
        Initial estimate of variance of eps(1) = .04
        Initial estimate of variance of eps(2) = .12
      Each  initial  estimate  may  optionally be coded with one of the
      forms:
         init options      (init options)        (options init)

      With NONMEM 7.3 (value,value...)xn is permitted, so that repeated
      inputs  of  $SIGMA  may  be entered easily.  Any initial value or
      group of initial values may be enclosed in parentheses  and  fol-
      lowed  by "xn", which means to replicate the values within paren-
      theses n times ("repeated value").

      The following options apply only to  a  single  initial  estimate
      (i.e.,  a  single 1x1 block) and must follow the initial estimate
      unless within parentheses.

      Option FIXED indicates that the variance is to be constrained  to
      be fixed to the given initial estimate.  (When FIXED appears any-
      where, then the block is described by  NONMEM  as  consisting  of
      separate blocks, each of dimension one.)

      Option  UNINT  is used with NONMEM 7.5.  UNINT is used during the
      Optimal Design Step to identify an eps as  uninteresting.   UNINT
      may be used anywhere that FIXED may be used.

      Option VARIANCE indicates that the initial estimate is understood
      to be a variance of the eps.  This is the default.

      Option STANDARD indicates that the initial estimate is understood
      to be a standard deviation of the eps.  May also be coded SD.

      An  initial  estimate  may  be 0 only if the variance or standard
      deviation is fixed to this estimate.

 2.   $SIGMA BLOCK(n) [ v11 v21 v22 v31 v32 v33 ... vn1 vn2 ... vnn ]

      This gives the initial estimates of all the elements of a  nondi-
      agonal ("full") block of SIGMA.  E.g.,
        $SIGMA BLOCK(2) .04 .002 .12
        Initial estimate of variance of eps(1) = .04
        Initial estimate of covariance of eps(2), eps(1) = .002
        Initial estimate of variance of eps(2) = .12

      Any  initial  value or group of initial values may be enclosed in |
      parentheses and followed by "xn", which means  to  replicate  the |
      values within parentheses n times ("repeated value").

      The  following  options  apply to the entire block and may appear
      anywhere among the list of initial estimates:

      FIXED indicates that the entire block is constrained to be  fixed
      to its initial estimate.

      Option  UNINT  is used with NONMEM 7.5.  UNINT is used during the
      Optimal Design Step to identify an eps as  uninteresting.   UNINT
      may be used anywhere that FIXED may be used.

      VARIANCE  indicates that all initial estimates given for diagonal
      elements are understood to be initial estimates of  variances  of
      eps.  This is the default.

      STANDARD  indicates that all initial estimates given for diagonal
      elements are understood to be initial estimates of standard devi-
      ations of epsilons.  May also be coded SD.

      COVARIANCE  indicates  that  all  initial estmates given for off-
      diagonal elements are  understood  to  be  initial  estimates  of
      covariances of epsilons.  This is the default.

      CORRELATON  indicates  that  all  initial estmates given for off-
      diagonal elements are understood to be initial estimates of  cor-
      relations of epsilons.

      CHOLESKY  indicates  that  the block is specified in its Cholesky
      form.

      Options VARIANCE or STANDARD may be combined with  COVARIANCE  or
      CORRELATON.

      Note  that  NONMEM converts all initial estimates to variance and
      covariances.  The values desplayed in the NONMEM  report  and  in
      the  raw  and  additional  output  files are always variances and
      covariances.

      Examples:

      The following describe the same block (within rounding errors):
      $SIGMA BLOCK(2) ; or $SIGMA VARIANCE COVARIANCE BLOCK(2)
      0.64
      -0.24 0.58
      $SIGMA STANDARD BLOCK(2)
      0.8
      -0.24 0.762
      $SIGMA STANDARD CORRELATION BLOCK(2)
      0.8
      -0.394 0.762
      $SIGMA VARIANCE CORRELATION BLOCK(2)
      0.64
      -0.394 0.58
      $SIGMA CHOLESKY BLOCK(2)
      0.8
      -0.3 0.7

      The (entire) initial estimate of the block must be positive defi-
      nite.   The only exception is when the entire initial estimate of
      the block is 0, in which case it must be fixed to this  estimate.
      Initial  estimates of some of the elements of the block may be 0,
      while initial estimates of some other elements  may  be  nonzero,
      but only in the case where the block is constrained to be of band
      symmetric form.  That is, given the diagonal and a group of  con-
      tiguous  subdiagonals symmetrically ocurring across the diagonal,
      the elements off both the diagonal and the subdiagonals are  con-
      strained  to be zero.  To specify the initial estimates of such a
      block, the initial estimates of those elements  that  are  to  be
      constrained  to  0  should be given as 0, while all other initial
      estimates should be given as nonzero.  E.g.,  with  these  struc-
      tures for $SIGMA BLOCK(3), the 0's are preserved:
       x
       0x
       00x

       x
       xx
       0xx

      With  NONMEM  7.3 if the initial estimate of a block is not posi- |
      tive definite because of rounding errors, a value will  be  added |
      to  the diagonal elements to make it positive definite. A message |
      in the NONMEM report file  will  indicate  that  this  was  done. |
      E.g.,                                                             |
      DIAGONAL SHIFT OF  1.1000E-03 WAS IMPOSED TO ENSURE POSITIVE DEF- |
      INITENESS

 3.   $SIGMA BLOCK(n) SAME(m)

      This describes a block whose initial estimates, as well as  final
      estimates,  are constrained to be equal to those of the preceding
      block. Values may not be given. "(n)" may be omitted.
      With NONMEM 7.3 (m) is permitted.  If (m) is present,  then  this
      record is equivalent to m identical records without (m).  E.g.,
      $SIGMA BLOCK(2) SAME(3)
      is equivalent to
      $SIGMA BLOCK(2) SAME
      $SIGMA BLOCK(2) SAME
      $SIGMA BLOCK(2) SAME

 4.
      $SIGMA BLOCK(n) VALUES(diag,odiag)

      This  supplies  initial  values for a block such that the initial
      estimates of the diagonal elements are all the same, specified by
      "diag",  and  the  initial estimates of the off-diagonal elements
      are all the same, specified by "odiag".  If present, VALUES  must
      follow  BLOCK.   Other  options  (such  as FIXED, CHOLESKY, VARI-
      ANCE,STANDARD,COVARIANCE,CORRELATON, UNINT) may follow VALUES  or
      be placed between BLOCK and VALUES.
      E.g.,
      $SIGMA BLOCK(6) VALUES(0.1,0.01)
      is the same as
      $SIGMA BLOCK(6)
      0.1
      0.01 0.1
      (0.01)x2 0.1
      (0.01)x3 0.1
      (0.01)x4 0.1
      (0.01)x5 0.1

      For fixed block (such as for sigma priors):
      $SIGMA BLOCK(6) FIX VALUES(0.15,0.0)

 5.
      $SIGMA label=value (NM75)

      The  symbolic  label substitution feature is new with NONMEM 7.5.
      This is a compact method of defining an EPS (an element of SIGMA)
      specifying  its  initial estimate, and specifying a label for the
      subscript for this element of SIGMA.  The label may be used as  a
      subscript  for  EPS  in  abbreviated code, and will also identify
      this element of SIGMA  in  the  NONMEM  output.   If  new  $SIGMA
      records  change the ordering, the abbreviated code does  not have
      to be changed.  For example, suppose the first element  of  SIGMA
      that is defined happens to be
      $SIGMA RSW=.06
      The NONMEM report will describe the relationship, e.g.,
      LABELS FOR EPS
      EPS(1)=EPS(RSW)
      and  EPS(RSW)  rather than EPS1 will appear in the NONMEM report.
      The abbreviated code can use this symbolic subscript  instead  of
      the numeric subscript.

      As  with  $OMEGA, $SIGMA and $THETA records (if elements of THETA
      are used) must be placed ahead of any records that use  the  sym-
      bolic label.

      Another example defines symbolic labels for a block of SIGMA:

      $SIGMA BLOCK(2)
      RSW=  0.3
      EX=   0.01 0.35

      Or, for diagonals,

      $SIGMA
      RSW= 0.3
      EX= 0.35

 6.
      $SIGMA BLOCK(n) NAMES (label1,...,labeln) VALUES (odiag,diag)  (NM75)
      With NONMEM 7.5, Symbolic label substitution may be specified for
      an entire block using the NAMES option.  This is a compact way of
      defining one or more epsilons with labels and, when combined with
      VALUES, with initial values. For example

      $SIGMA BLOCK(2) NAMES(RSW,EX) VALUES(0.03,0.01)

      This is equivalent to

      $SIGMA BLOCK(2)
      RSW=  0.03
      EX=   0.01 0.03

 If both are present, VALUES() must come after NAMES().

 REFERENCES: Guide I, section C.3.4.6 , D.5.2 , D.5.3 
 REFERENCES: Guide IV, section III.B.10 , V.C.6 
 REFERENCES: Guide V, section 9.3 
