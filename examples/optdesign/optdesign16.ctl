
$PROB RUN# Example 1 (from samp5l)
$INPUT C SET ID JID TIME  DV=CONC AMT=DOSE RATE EVID MDV CMT TSTRAT TMIN TMAX
$DATA optdesign2.csv IGNORE=C

$SUBROUTINES ADVAN3 TRANS4

$PK
; The thetas are MU modeled.  
; Best that there is a linear relationship between THETAs and Mus
; The linear MU modeling of THETAS allows them to be efficiently 
; Gibbs sampled.

MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
CL=DEXP(MU_1+ETA(1))
V1=DEXP(MU_2+ETA(2))
Q=DEXP(MU_3+ETA(3))
V2=DEXP(MU_4+ETA(4))
S1=V1
CLN=THETA(1)
V1N=THETA(2)
QN=THETA(3)
V2N=THETA(4)

$ERROR
Y = F + F*EPS(1)+EPS(2)

; Initial values of THETA
$THETA 
 1.68338E+00  1.58812E+00  8.12710E-01  2.37436E+00 


;INITIAL values of OMEGA
;$OMEGA BLOCK(4) VALUES(0.0225,0.001)
$OMEGA (0.0225 FIXED)X4

;Initial value of SIGMA
$SIGMA 
( 0.0225)
( 0.0001 FIXED)

$PRIOR NWPRI PLEV=0.9999

$THETAP  (1.68338E+00 FIXED) (1.58812E+00 FIXED) (8.12710E-01 FIXED) (2.37436E+00 FIXED)

$THETAPV BLOCK(4) FIXED VALUES(0.03,0.001)

$DESIGN STGR EOPTD=1 FIMTYPE=1 SEED=1222345 DESEL=TIME DESELSTRAT=TSTRAT DESELMIN=TMIN DESELMAX=TMAX 
           MAXEVAL=9999 SIGL=12 nohabort PRINT=10 NOPRIOR=1
$TABLE ID CLN V1N QN V2N TIME EVID MDV DV NOPRINT NOAPPEND FILE=optdesign16.tab  FORMAT=S1PE23.16
