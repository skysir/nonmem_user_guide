$PROB RUN# example6 (from r2compl)
$INPUT C SET ID JID TIME DV=CONC DOSE=AMT RATE EVID MDV CMT WT
$DATA finetest5i.csv IGNORE=C

$FINEDATA tstart=FIRST TSTOP=LAST TDELTA=0.5 AXIS=TIME(LIN) CMT=1,3
          file=finetest5.csv