$PROB RUN# Example 1 (from samp5l)
$INPUT C SET ID JID TIME DATE=DROP DV=CONC AMT=DOSE RATE EVID MDV CMT CLX V1X QX V2X SDIX SDSX
$DATA example1td.csv IGNORE=C LAST20=80

$FINEDATA tstart=0 TSTOP=100 NEVAL=250 AXIS=TIME(LIN) CMT=1
          file=example1tdf.csv