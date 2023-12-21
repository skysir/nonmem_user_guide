$PROB RUN# Example 1 (from samp5l)
$INPUT ID TIME  DV=CONC AMT=DOSE RATE EVID MDV CMT
$DATA tdist11.csv IGNORE=C

$SUBROUTINES ADVAN3 TRANS4

$PK
MU_1=THETA(1)
MU_2=THETA(2)
MU_3=THETA(3)
MU_4=THETA(4)
NU=4.0
CHISQ=SQRT( (ETA(5)*ETA(5)+ETA(6)*ETA(6)+ETA(7)*ETA(7)+ETA(8)*ETA(8))/NU )
CL=EXP(MU_1+ETA(1)/CHISQ)
V1=EXP(MU_2+ETA(2)/CHISQ)
Q= EXP(MU_3+ETA(3)/CHISQ)
V2=EXP(MU_4+ETA(4)/CHISQ)
S1=V1

$ERROR
Y = F + F*EPS(1)

$THETA 1.68338E+00  1.58811E+00  8.12694E-01  2.37435E+00  
$OMEGA BLOCK(4)
0.03
0.01  0.03 
-0.006 0.01  0.03 
0.01 -0.006  0.01 0.03

$OMEGA (1.0 FIXED) (1.0 FIXED) (1.0 FIXED) (1.0 FIXED)

$SIGMA 
0.01

;$EST METHOD=ITS INTERACTION MAXEVAL=9999 PRINT=5 NOHABORT SIGL=9 CTYPE=3 NITER=200 NONINFETA=1 MCETA=10
$EST METHOD=IMP INTERACTION MAXEVAL=9999 PRINT=1 NOHABORT ISAMPLE=1000 NITER=200 SIGL=9 DF=2 RANMETHOD=3S2P
     CTYPE=3 MCETA=10
$EST METHOD=1 INTERACTION MAXEVAL=9999 PRINT=1 NOHABORT NSIG=3 SIGL=9 NONINFETA=1 SLOW MCETA=10
$COV MATRIX=R UNCONDITIONAL
