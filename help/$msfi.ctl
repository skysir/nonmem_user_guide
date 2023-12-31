


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $MSFI                                |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Gives the name of an input Model Specification File
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $MSFI   filename [NORESCALE|RESCALE] [NPOPETAS[=n]]
         [ONLYREAD]  [NOMSFTEST|MSFTEST]
         [VERSION=n]                                                    |
         [NEW]                                                          |

 SAMPLE:
 $MSFI   MSF13

 DISCUSSION:
 This record gives the name of an input Model Specification File.

 OPTIONS:

 filename
      Name  of the Model Specification File.  Must be the first option.
      Filename may not contain embedded spaces.  If it contains commas,
      semicolons,  or parentheses, then it must be surrounded by single
      quotes ' or double quotes ".  Filename  may  also  contain  equal
      signs  if  it  is enclosed in quotes.  If filename is the same as
      any option of the $MSFI record, it must be  enclosed  in  quotes.
      Filename may contain at most 71 characters.

 NORESCALE
      If  a  search is continued, it is continued just as it would have
      been in the previous run.  This option has no effect if the Esti-
      mation Step is not implemented or if the step is implemented, but
      using different options on the $ESTIMATION record from those used
      with the previous run.  This is the default option.

 RESCALE
      The  search is to be restarted with initial estimates taken to be
      the final estimates from the previous run, so that a UCP value of
      0.1  now  corresponds  to a final estimate from the previous run.
      This option may be used only if the search with the previous  run
      terminated  successfully and all the options used on the $ESTIMA-
      TION record coincide  with  those  used  with  the  previous  run
      (except the MAXEVAL option).

 NPOPETAS[=n]
      If  the  data  are  to  be  understood to be population data, the
      option NPOPETAS may be needed, and then the part "=n",  n>0,  may
      also  need  to  be  coded.  The integer n should be the number of
      population eta variables.  The part "=n"  is  required  when  the
      $ERROR  block  describes a model that changes the value of F, and
      there is no $PK block.  It is also required when  variables  A(n)
      are  used  in  the  $PK  abbreviated code.  If the data are to be
      understood to be single-subject data, the option  NPOPETAS=0  may
      be needed.
      May also be coded POPETAS.

 ONLYREAD
      This option should be used if and only if the model specification
      file is being used to convey prior information  to  a  succeeding
      problem  in  the control stream (See TNPRI).  No task records may
      appear in a problem specification containing a $MSFI record  with
      this option.

 When  the  $MSFI  record  is  used in a problem specification, $THETA,
 $OMEGA and $SIGMA records should not appear in that specification.  It
 is  permissible  for  the  data  file to differ from the original data
 file, but usually this is not done.

 MSFTEST
      Prevents an MSF file from being utilized in a subsequent  control
      stream  file  or  problem  if  there  were  errors.   This is the
      default.

 NOMSFTEST
      Sometimes the MSFI error check is too strict.  This  occurs  par-
      ticularly  when  using classical NONMEM methods.  To turn of MSFI
      error checking, use NOMSFTEST.

 VERSION[=n]
      Read in MSF files from  previous  versions.  With  NONMEM  7.5.0, |
      choices are:                                                      |
      $MSFI myfile.msf VERSION=7.4.0                                    |
      $MSFI myfile.msf VERSION=7.3.0                                    |
      $MSFI myfile.msf VERSION=7.2.0                                    |
      $MSFI myfile.msf VERSION=7.1.2                                    |
      $MSFI myfile.msf VERSION=7.1.0                                    |
      $MSFI myfile.msf VERSION=6.2                                      |
      $MSFI myfile.msf VERSION=6.1                                      |

 NEW  When  the problem that created the MSF file has successfully com- |
      pleted, calling for a resumed or new estimation is prevented when |
      the  method is FO/FOCE/Laplace. To allow analysis to continue, or |
      to allow an analysis on a new data set, resuming from  the  final |
      parameters  of the MSF file, use the option NEW: $MSFI myfilename |
      NEW.

 With NONMEM 7.3 and later, when MSF or MSFO option is used to  specify
 an MSFO file in the $EST record e.g.,
 $EST ... MSFO=msfroot.msf
 then  in  addition  to  the  main  MSF  file  msfroot.msf,  file  msf-
 root_ETAS.msf containing individual etas will also  be  produced,  and
 provide additional information when a $MSFI record is used in a subse-
 quent problem or control stream.  This is referred to  as  an  "extra"
 msf  file.   If  the  Covariance  Step is also implemented, files msf-
 root_RMAT.msf and msfroot_SMAT.msf containing intermediate information
 on  the R matrix and S matrix will also be produced.  These files pro-
 vide information when a $MSFI record along  with  a  $COV  ...  RESUME
 record is used in a subsequent problem or control stream.

 The use of an extension, e.g., .msf, is optional.  If the _ETA file is
 not present, NONMEM issues a warning:
 WARNING: EXTRA MSF FILE COULD NOT BE OPENED: c5msf2x_ETAS

 There is no warning if _SMAT and/or _RMAT are not present.

 REFERENCES: Guide I, section B.3 
 REFERENCES: Guide II, section E.3 
 REFERENCES: Guide III, section 3.2 
 REFERENCES: Guide IV, section III.B.12 
 REFERENCES: Guide Introduction_7
