


 +--------------------------------------------------------------------+
 |                                                                    |
 |                               $SIZES                               |
 |                                                                    |
 +--------------------------------------------------------------------+

 MEANING: Array sizes for NONMEM and PREDPP
 CONTEXT: NM-TRAN Control Record

 USAGE:
 $SIZES  [constant=value] [constant=value] ...

 SAMPLE:
 $SIZES LIM1=30000 MAXFCN=2000000 NO=500

 DISCUSSION:
 $SIZES is optional.  If present, it must precede the first $PROBLEM or
 $SUPER record.

 Certain constants are used in NM-TRAN, NONMEM and PREDPP.  With NONMEM
 7.2  and  higher, the user may override many of the constants with the
 $SIZES record.

 See the discusion of sizes for a discussion of how these constants are
 determined, and how they are communicated to NONMEM and PREDPP.
 (See SIZES FSIZES prsizes)

 Any  non-zero  value specified on the $SIZES record overrides both the
 default and the value that NM-TRAN would have specified.  (A value  of
 0  is  ignored.)   As  of  NONMEM  7.3, as an alternative to modifying |
 sizes.f90 to very large maximum sizes, you can tell NMTRAN the maximum |
 size  that may be needed by specifying a $SIZES constant as a negative |
 value. Thus, a user can give NMTRAN permission to deal with all  prob- |
 lems  that  have data input files that have up to 1000 data items, and |
 up to 150 etas and epsilons, and up to 200 thetas, by the following:   |

 $SIZES PD=-1000 LVR=-150 LTH=-200                                      |

 but the values of these constants when the NONMEM executable  is  con- |
 structed  will  be  only what is needed for the particular problem. In |
 contrast,                                                              |

 $SIZES PD=1000 LVR=150 LTH=200                                         |

 will result in sizing the NONMEM executable  with  these  values,  and |
 won't  make  a  "tailor  fit".  This would result in a very large exe- |
 cutable regardless of the model size. Thus,                            |
 $SIZES PD=-1000                                                        |
 tells NMTRAN that you may need as many as 1000 data items  in  a  data |
 file, whereas                                                          |
 $SIZES PD=1000                                                         |
 tells NMTRAN that you need exactly that size.

 List of $SIZES Record Options and Their Default Values

     LTH=100
     LVR=30
     LVR2=20
     NO=250
     MMX=10
     LNP4=4000
     LSUPP=4050
     LIM7=2
     LWS3=9000
     MAXIDS=10000
     LIM1=10000
     LIM2=100000
     LIM3=10000
     LIM4=1000
     LIM5=200
     LIM6=400
     LIM8=200
     LIM11=25
     LIM13=1000
     LIM15=1000
     LIM16=400
     MAXRECID=200
     PC=30
     PCT=30
     PIR=700
     PD=50
     PDT=50
     PAL=50
     MAXFCN=1000000
     DIMTMP=500
     DIMCNS=500
     DIMNEW=1000
     DIMVRB=200
     PL=10
     NPOPMIXMAX=10
     MAXOMEG=70
     MAXPTHETA=90
     MAXITER=210
     ISAMPLEMAX=10
     MAXSIDL=0
     PNM_MAXNODES=100
     MAXNRDS=PC
     PAST_SIZE=4000

 Additional constants that may be set with $SIZES:

 LADD_MAX (See resource/SIZES.f90).

 (See SIZES FSIZES prsizes)

 REFERENCES: Guide Introduction_7
