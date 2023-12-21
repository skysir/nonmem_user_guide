


 +--------------------------------------------------------------------+
 |                                                                    |
 |                    CTLO CTUP PROBABILITY: PR_CT                    |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: NONMEM read-only global variables
 CONTEXT: User-supplied routines

 USAGE:
      USE ROCM_REAL, ONLY PR_CT

 GLOBAL DECLARATION:
      USE SIZES, ONLY: DPSIZE
      REAL(KIND=DPSIZE) :: PR_CT

 DISCUSSION:

 When  with  a given data record, either of the limits CTLO or CTUP are
 set, thus defining an interval of values  comprising  one  of  several
 categories  equated  with  the possible values of a potential observa-
 tion, then during a copying pass (and during ICALL=5 and 6), PR_CT  is
 the  estimated probability that an observation will be of the category
 in question.
 (See CTLO, CTUP)
 (See copying_block, expectation, data_average)

 If the mean and variance of the intra-individual model for a potential
 observation  are  specified, and if the limits CTLO or CTUP are set, a
 value of PR_CT will be computed,  whether  the  record  has  MDV=0  or
 MDV=1.  If neither CTLO nor CTUP are set, the value PR_CT will be 1.

 PR_CT  may  be  used  as a right-hand quantity in abbreviated code for
 $PRED, $PK, $ERROR, and $INFN blocks.

 Location prior to NONMEM 7: rocm45

 REFERENCES: None.
