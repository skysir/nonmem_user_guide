$PROB RUN# example6 (from r2compl)
$INPUT C SET ID JID TIME DV=CONC DOSE=AMT RATE EVID MDV CMT
$DATA example6.csv IGNORE=C

$FINEDATA TSTART=0.01 TSTOP=50 NEVAL=100 AXIS=TIME(LOG) CMT=1,3
          FILE=example6e.csv